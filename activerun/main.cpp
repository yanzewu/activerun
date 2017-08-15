#include "activerun.h"


int main(int argc, char* argv[]) {

    if (argc < 2) {
        printf("Usage activerun [inputfile]");
        return 1;
    }

    // prepare system

    InputParameter input;
	if (input.read_input(argv[1])) {
		return 1;
	}


    DataFile datafile;
    const char* data_filename = "data.gel";
    datafile.read_data(data_filename);

    System system;
    system.init(datafile);

    State state_init;
    state_init.init(datafile);

    Context context;
    context.init_force_buffer(system, 3);
	context.timestep = input.dt;

    // prepare forces

    BrownianForce force_brown;
    force_brown.init(input, system);
    double brown_timescale = input.zeta / input.kT;
    
    SwimForce force_swim;
    bool has_swim = std::count(system.atom_type.begin(), system.atom_type.end(), 1) > 0;
    double swim_timescale = 100.0;
    if (has_swim) {
        force_swim.init(input, system);
		force_swim.update_cache(system, context);
    }

    MorseForce force_morse;
    force_morse.init(input, system);

//    context.timestep = std::min(brown_timescale, swim_timescale) * input.dt;
    context.total_steps = (size_t)(input.stepstR / context.timestep);
    double ave_atom_size = 1.0;
    context.init_neighlist(system.box, input.mincellL * std::max(force_morse.cutoff_global, ave_atom_size));
    context.init_pbc(system, true);
	context.init_multicore(input.np, 2);

	force_morse.init_mpi(context.thread_num[2]);

    // prepare integrator

    LangevinIntegrator integrator;
    integrator.compute_temperature = true;
    integrator.init(system, context);

    // prepare runtime

    State state = state_init;

    // prepare dumper

    TrajDumper trajdumper("outputconfig.lammpstrj");
    LineDumper swim_stress_dumper("swimstress.output", {"stressPa"});
    LineDumper total_stress_dumper("totalstress.output", {});

	ThermoStat thermostat;
	thermostat.init(context, state_init);
	thermostat.compute_pressure = { 1, 1, 1 };
	thermostat.compute_temperature = false;

    // run
    integrator.update_cache(context);
	thermostat.update_cache(system);
    force_brown.update_cache(system, context);
    force_morse.update_cache(system, context);

    trajdumper.dump(system, state, 0);

    std::vector<double> pressure(3, 0.0);
    double temperature;

	FixedThreadPool pool;
	pool.init(4);

    for (context.current_step = 0; context.current_step < context.total_steps; context.current_step++) {

//		integrator.update_first_half(state);
        try {
            context.pbc.update(state.pos);
        }
        catch (const std::out_of_range&) {
            fprintf(stderr, "At step %d\n", context.current_step);
			dump(state, context);
			trajdumper.dump(system, state, context.current_step);
            return 1;
        }
        context.neigh_list->build_from_pos(state.pos);
		for (auto& fb : context.force_buffer) {
			memset(&fb[0], 0, sizeof(Vec) * fb.size());
		}

	//	pool.submit(std::bind(&BrownianForce::update, &force_brown, std::ref(state), std::ref(context.force_buffer[0])));

        if (has_swim && context.current_step >= input.swimstart) {
	//		pool.submit(std::bind(&SwimForce::update, &force_swim, std::ref(state), std::ref(context.force_buffer[1])));
            force_swim.update(state, context.force_buffer[1]);
        }
		force_morse.update(pool, state, context.pbc, *context.neigh_list);
		pool.start_all();
        force_brown.update(state, context.force_buffer[0]);
		pool.wait();
		force_morse.update_later(context.force_buffer[2]);

        integrator.update(state, context);
		thermostat.temperature_cache = integrator.update_temperature();
        if ((context.current_step + 1) % input.output == 0) {
			context.pbc.update_image(state.pos);
			thermostat.update(context, integrator.velocity_cache);
			printf("%zd\t% .5e\t% .5e\t% .5e\t% .3f\n", 
				context.current_step + 1,
				thermostat.pressure_cache[0], 
				thermostat.pressure_cache[1], 
				thermostat.pressure_cache[2], 
				thermostat.temperature_cache);
            trajdumper.dump(system, state, context.current_step + 1);
 //           swim_stress_dumper.dump(state, step);
  //          total_stress_dumper.dump(state, step);
        }
    }
}