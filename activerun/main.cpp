#include "activerun.h"
#include "timer.h"
#include "Serializer.h"
#include "old_input.h"

void set_default_value(Serializer& config) {
    
    // global
    config["data_file"] = std::string("data.gel");
    config["is_restart"] = false;
    config["global_seed"] = 0;
    config["current_step"] = 0;

    // run
    config["run"] = json({
        {"step", 400000},
        {"timestep", 5e-6}
    });

    // acceleration
    config["util"] = json({
        {"cell_size", 2.0},
        {"np", 4},
        {"neighlist_step", 1}
    });

    // dump
    config["dump"] = json({
        {"dump_file", "outputconfig.lammpstrj"},
        {"dump_step", 10000}
    });

    // thermo
    config["thermo"] = json({
        {"using_thermo", true},
        {"thermo_file", "pressure.output"},
        {"thermo_step", 10000},
        {"thermo_start", 0},
        {"average_step", 1},
        {"average_range", 0},
        {"compute_temp", true}
    });

    // time
    config["timer"] = json({
        {"sample_step", 100}
    });

    // restart
    config["restart"] = json({
        {"write_restart", true},
        {"restart_file", "restart.rst"},
        {"restart_data", "restart.data"}
    });

    // forces
    config["BrownianForce"] = json({ { "temp", 1.0 },{ "neta", 1.0 } });
    config["SwimForce"] = json({ 
        { "swim_start", 0},
        { "swim_temp", 1.0 },
        { "neta", 1.0 }, 
        { "brownian", true},
        { "PeR", 0.01},
        { "init_seed", 0}
    });
    config["MorseForce"] = json({
        { "kappa", 30.0 },
        { "Um", 5.0 },
        { "cutoff", 1.25}
    });
}

int read_input(int argc, char* argv[], Serializer& config) {

    // prepare system
    if (argc < 2) {
        printf("Usage: activerun [inputfile] (datafile)\n");
        printf("activerun -r [restart]\n");
        printf("activerun -i [input json]\n");
        return 1;
    }

    if (strcmp(argv[1], "-r") == 0) {
        RestartFile restart;
        restart.read_restart(argv[2]);
        if (strncmp(restart.input_name + strlen(restart.input_name) - 4, "json", 4) == 0) {
            config.read_file(restart.input_name);
            config["input_file"] = std::string(restart.input_name);
        }
        else {
            return read_legacy_input(restart.input_name, config);
        }

        config["is_restart"] = true;
        config["data_file"] = std::string(restart.data_name);
        config["current_step"] = restart.current_step;

        config["dump"]["dump_file"] = std::string(restart.output_name);
        if (strlen(restart.thermo_name) > 0) {
            config["thermo"]["using_thermo"] = true;
            config["thermo"]["thermo_file"] = std::string(restart.thermo_name);
        }
        else {
            config["thermo"]["using_thermo"] = false;
        }
    }
    else if (strcmp(argv[1], "-i") == 0) {
        config.read_file(argv[2]);
        config["input_file"] = std::string(argv[2]);
    }
    else {
        if (argc >= 3) {
            config["data_file"] = std::string(argv[2]);
        }
        return read_legacy_input(argv[1], config);
    }

    return 0;
}

int main(int argc, char* argv[]) {

    Serializer config;

    set_default_value(config);

    if (read_input(argc, argv, config)) {
        fprintf(stderr, "Error reading input\n");
        return 1;
    }


    json util_param = config["util"];
    json dump_param = config["dump"];
    json thermo_param = config["thermo"];
    json restart_param = config["restart"];
    json swim_param = config["SwimForce"];
    json timer_param = config["timer"];

    Dict fix_brownian_param = config.get_dict("BrownianForce");
    Dict fix_swim_param = config.get_dict("SwimForce");
    Dict pair_param = config.get_dict("MorseForce");
    Dict integrator_param = { { "compute_temp", config["thermo"]["using_thermo"] && config["thermo"]["compute_temp"] ? 1.0 : 0.0 } };
    Dict time_param = config.get_dict("run");

    bool is_restart = config["is_restart"];
    std::string data_file = config["data_file"];

    DataFile datafile;
    datafile.read_data(data_file.c_str());

    System system;
    system.init(datafile);

    State state_init;
    state_init.init(datafile);

    Context context;

    // prepare forces

    context.init_force_buffer(system, 3); // 3 forces

    BrownianForce force_brown;
	force_brown.init(config.get_dict("BrownianForce"), system);
    
#ifndef THREE_DIMENSION
    SwimForce force_swim;
#else
    SwimForce3d force_swim;
#endif // !THREE_DIMENSION

    bool has_swim = std::count(system.atom_type.begin(), system.atom_type.end(), 1) > 0;
    size_t swim_start = swim_param["swim_start"];

	if (has_swim) {
        force_swim.init(fix_swim_param, system);
	}

    MorseForce force_morse;
	force_morse.init(pair_param, system);

    double cell_size = util_param["cell_size"];
    int np = util_param["np"];
    int neighlist_step = util_param["neighlist_step"];

	context.init_timestep(time_param);
    context.init_pbc(system, true);
	double max_atom_size = *std::max_element(system.get_attr("size").begin(), system.get_attr("size").end());
    context.init_neighlist(system.box, cell_size * std::max(force_morse.max_cutoff(), max_atom_size));
	context.init_multicore(np, 2); // force 2 is paired force

	// mpi

	if (has_swim && np > 2) {
		context.thread_num[2]--;
		context.thread_num[0] = 1;
	}

	force_brown.init_mpi(context.thread_num[0]);
	if (has_swim) {
		force_swim.init_mpi(context.thread_num[1]);
	}
	force_morse.init_mpi(context.thread_num[2]);

    // integrator

    LangevinIntegrator integrator;
    integrator.init(integrator_param, system, context);

    // dump

    std::string dump_file = dump_param["dump_file"].get<std::string>().c_str();
    size_t dump_step = dump_param["dump_step"].get<size_t>();

    TrajDumper trajdumper(dump_file.c_str(), is_restart);

    // thermo

    size_t step_begin = config["current_step"];

    std::string thermo_file = thermo_param["thermo_file"].get<std::string>().c_str();
    bool using_thermo = thermo_param["using_thermo"];
    size_t thermo_step = thermo_param["thermo_step"];
    size_t thermo_start = thermo_param["thermo_start"];
    int thermo_ave_range = thermo_param["average_range"];
    int thermo_ave_step = thermo_param["average_step"];

    if (thermo_start < step_begin) {
        printf("Warning: thermo_start too small (%zd), using first frame (%zd) instead.\n", thermo_start, context.current_step);
        thermo_start = step_begin;
    }

    std::vector<double> pressure(3, 0.0);
    double temperature;
    ThermoStat thermostat;
    int count = 0;

    LineDumper thermodumper(thermo_file.c_str(), { "P_Brown", "P_Swim", "P_Morse", "Temp" }, true, is_restart);

    // cache

    printf("\nInitializing...\n");

    integrator.update_cache(system, context);
    if (using_thermo) {
	    thermostat.update_cache(system);
    }
    force_brown.update_cache(system, context);
    if (has_swim) {
		force_swim.update_cache(system, context);
    }
    force_morse.update_cache(system, context);

    // run

	printf("\nStart running session...\n\n");

    State state = state_init;

    if (!is_restart) {

        trajdumper.dump(system, state, 0);

        if (using_thermo) {
	        thermodumper.dump_head();
        }
    }

    // random seed

    int global_seed = config["global_seed"];
	set_random_seed(global_seed);


    // timer

    Timer timer;
    timer.init({ "NeighList", "Force", "Integrator" });
    int timer_step = timer_param["sample_step"];

    for (context.current_step = step_begin; context.current_step < context.total_steps; context.current_step++) {

        if (using_thermo) {
            if (context.current_step == thermo_start) {
                context.pbc.reset_location();
                context.pbc.update_image(state.pos);

                thermostat.init(context, state);
                thermostat.compute_pressure = { 1, 1, 1 };
                thermostat.compute_temperature = false;
                integrator.compute_temperature = true;
                integrator.update_cache(system, context);
                temperature = 0.0;
            }
        }

        if (context.current_step % timer_step == 0) timer.set_checkpoint(-1);

		if (context.current_step % neighlist_step == 0) {
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
		}

        if (context.current_step % timer_step == 0) timer.set_checkpoint(0);

		for (auto& fb : context.force_buffer) {
			memset(&fb[0], 0, sizeof(Vec) * fb.size());
		}
        //try {
		force_morse.mp_update(context.pool, state, context);    // morse force must update first
        /*
        }
        catch (const std::runtime_error&) {
            trajdumper.dump(system, state, context.current_step / dump_step * dump_step);
            dump_snapshot(state, context);
            return 1;
        }*/
        force_brown.mp_update(context.pool, state, context.force_buffer[0]);
        if (has_swim && context.current_step >= swim_start) {
			force_swim.mp_update(context.pool, state, context.force_buffer[1]);
        }
		context.pool.wait();
		force_morse.update_later(context.force_buffer[2]);

        if (context.current_step % timer_step == 0) timer.set_checkpoint(1);

        integrator.update(state, context);

        if (context.current_step % timer_step == 0) timer.set_checkpoint(2);

        if ((context.current_step + 1) % dump_step == 0) {
            trajdumper.dump(system, state, context.current_step + 1);
        }

        if (using_thermo) {
            if (context.current_step + 1 >= thermo_start + thermo_ave_range / 2) {
                int step_overhead = (context.current_step + 1 + thermo_ave_range / 2 - thermo_start) % thermo_step;
                if (step_overhead <= thermo_ave_range) {
                    if (step_overhead % thermo_ave_step == 0) {
                        count++;
                        context.pbc.update_image(state.pos);
                        thermostat.update(context, integrator.velocity_cache);
                        for (size_t i = 0; i < pressure.size(); i++) {
                            pressure[i] += thermostat.pressure_cache[i];
                        }
                        temperature += integrator.update_temperature();
                    }
                }
                
                if (step_overhead == thermo_ave_range) {
                    thermodumper.dump({
                        pressure[0] / count,
                        pressure[1] / count,
                        pressure[2] / count,
                        temperature / count },
                        context.current_step + 1 - thermo_ave_range/2);
                    memset(&pressure[0], 0, pressure.size() * sizeof(double));
                    temperature = 0.0;
                    count = 0;
                }
            }
        }
    }

    timer.print((context.current_step - step_begin) / timer_step, context.current_step - step_begin);

    bool write_restart = restart_param["write_restart"];

	if (write_restart) {
        RestartFile restart;

        std::string input_file = config["input_file"];
        std::string restart_file = restart_param["restart_file"];
        std::string restart_data = restart_param["restart_data"];

		state.write_data(datafile);
		datafile.write_data(restart_data.c_str());

        restart.current_step = context.current_step;
        strcpy(restart.data_name, restart_data.c_str());
        strcpy(restart.input_name, input_file.c_str());
        strcpy(restart.output_name, dump_file.c_str());
        if (using_thermo) {
            strcpy(restart.thermo_name, thermo_file.c_str());
        }
        restart.write_restart(restart_file.c_str());
		printf("\nRestart written to %s\n", restart_file.c_str());
	}

	return 0;
}