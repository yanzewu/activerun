#include "activerun.h"


int main(int argc, char* argv[]) {


    // prepare system
    if (argc < 2) {
        printf("Usage activerun [inputfile] (datafile)");
        return 1;
    }

    InputParameter input;
	if (input.read_input(argv[1])) {
		return 1;
	}
    const char* data_filename = argc > 2 ? argv[2] : "data.gel" ;
	const char* dump_filename = "outputconfig.lammpstrj";
	char restart_filename[256];
	strcpy(restart_filename, data_filename);
	strcat(restart_filename, ".restart");

	bool write_restart = true;

	Dict fix_brownian_param = {
		{"temp", input.kT},
		{"neta", input.viscosity} };
	Dict fix_swim_param = {
		{"swim_temp", input.kT},
		{"neta", input.viscosity},
		{"brownian", input.brownrot ? 1.0 : 0.0},
		{"PeR", input.PeR},
		{"init_seed", 0.0}
	};
	Dict pair_param = {
		{ "kappa", input.kappa },
		{ "Um", input.Um },
		{ "cutoff", input.Rg * 2 + 1.0 } };
	Dict integrator_param = {
		{"compute_temp", 1.0}
	};
	Dict time_param = {
		{"timestep", input.dt},
		{"steps", input.stepstR / input.dt }
	};
	Dict mpi_param = {
		{"np", 4.0}
	};


    DataFile datafile;
    datafile.read_data(data_filename);

    System system;
    system.init(datafile);

    State state_init;
    state_init.init(datafile);

    Context context;

    // prepare forces

    context.init_force_buffer(system, 3);
	//std::vector<std::string> force_names = { "Brownian", "Swim", "Morse" };

    BrownianForce force_brown;
	force_brown.init(fix_brownian_param, system);
    
    SwimForce3d force_swim;
    bool has_swim = std::count(system.atom_type.begin(), system.atom_type.end(), 1) > 0;
	if (has_swim) {
        force_swim.init(fix_swim_param, system);
	}

    MorseForce force_morse;
	force_morse.init(pair_param, system);

	// integrator+

	context.init_timestep(time_param);
    context.init_pbc(system, true);
	double max_atom_size = *std::max_element(system.get_attr("size").begin(), system.get_attr("size").end());
    context.init_neighlist(system.box, input.mincellL * std::max(force_morse.cutoff_global, max_atom_size));

    LangevinIntegrator integrator;
    integrator.init(integrator_param, system, context);

	// mpi

	context.init_multicore((int)mpi_param.get("np", 1.0), 2);
	if (has_swim) {
//		context.thread_num[2]--;
	//	context.thread_num[0] = 1;
	}

	force_brown.init_mpi(context.thread_num[0]);
	if (has_swim) {
		force_swim.init_mpi(context.thread_num[1]);
	}
	force_morse.init_mpi(context.thread_num[2]);

    // prepare runtime

    State state = state_init;

    // prepare dumper


    std::vector<double> pressure(3, 0.0);
    double temperature;
	ThermoStat thermostat;
	thermostat.init(context, state_init);
	thermostat.compute_pressure = { 1, 1, 1 };
	thermostat.compute_temperature = false;

    // run

    integrator.update_cache(system, context);
	thermostat.update_cache(system);
    force_brown.update_cache(system, context);
    if (has_swim) {
		force_swim.update_cache(system, context);
    }
    force_morse.update_cache(system, context);

	printf("\nStart running session...\n\n");

    TrajDumper trajdumper(dump_filename);
    trajdumper.dump(system, state, 0);

    LineDumper thermodumper("totalstress.output", {"P_Brown", "P_Swim", "P_Morse", "Temp"}, true);
	thermodumper.dump_head();

	srand(0);
	int64_t time_neighlist = 0, time_force = 0, time_integrator = 0, time_total = 0, time_last = clock();

    for (context.current_step = 0; context.current_step < context.total_steps; context.current_step++) {

		time_last = clock();
//		if (context.current_step % 2 == 0) {
			try {
				context.pbc.update(state.pos);
			}
			catch (const std::out_of_range&) {
				fprintf(stderr, "At step %zd\n", context.current_step);
				dump_snapshot(state, context);
				trajdumper.dump(system, state, context.current_step);
				return 1;
			}
			context.neigh_list->build_from_pos(state.pos);

//		}

		time_neighlist += (time_total = clock()) - time_last;
		time_last = time_total;

		for (auto& fb : context.force_buffer) {
			memset(&fb[0], 0, sizeof(Vec) * fb.size());
		}
        try {
		force_morse.mp_update(context.pool, state, context);

        }
        catch (const std::runtime_error&) {
            trajdumper.dump(system, state, context.current_step / input.output * input.output);
            dump_snapshot(state, context);
            return 1;
        }
        force_brown.mp_update(context.pool, state, context.force_buffer[0]);
        if (has_swim && context.current_step >= input.swimstart) {
			force_swim.mp_update(context.pool, state, context.force_buffer[1]);
        }
	//	context.pool.start_all();
		context.pool.wait();
		force_morse.update_later(context.force_buffer[2]);

		time_force += (time_total = clock()) - time_last;
		time_last = time_total;

        integrator.update(state, context);

		time_integrator += (time_total = clock()) - time_last;
		time_last = time_total;

        if ((context.current_step + 1) % input.output == 0) {
			thermostat.temperature_cache = integrator.update_temperature();
            trajdumper.dump(system, state, context.current_step + 1);
			context.pbc.update_image(state.pos);
			thermostat.update(context, integrator.velocity_cache);
			thermodumper.dump({
				thermostat.pressure_cache[0],
				thermostat.pressure_cache[1],
				thermostat.pressure_cache[2],
				thermostat.temperature_cache },
				context.current_step + 1);
        }
    }

	time_total = clock();
	printf("\n\nTotal time: %.2fs\n", (double)time_total / CLOCKS_PER_SEC);
	printf("Force time:         %.2fs (%.2f%%)\n", (double)time_force / CLOCKS_PER_SEC, 100.0 * time_force / time_total);
	printf("Neighbourlist time: %.2fs (%.2f%%)\n", (double)time_neighlist / CLOCKS_PER_SEC, 100.0 * time_neighlist / time_total);
	printf("Integrator time:    %.2fs (%.2f%%)\n", (double)time_integrator / CLOCKS_PER_SEC, 100.0 * time_integrator / time_total);

	if (write_restart) {
		state.write_data(datafile);
		datafile.write_data(restart_filename);
		printf("\nRestart written to %s\n", restart_filename);
	}

	return 0;
}