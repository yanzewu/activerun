
/* Written by Yanze Wu @2017 */

#include "Compute.h"
#include "Context.h"
#include "Dumper.h"
#include "Force.h"
#include "Input.h"
#include "Integrator.h"
#include "Serializer.h"
#include "System.h"
#include "Thermo.h"
#include "timer.h"
#include "vecutil.h"

#ifdef READ_LEGACY_INPUT
#include "old_input.h"
#endif 

MultiDumper* logger;

void set_default_value(Serializer& config) {
    
    // global
    config["data_file"] = std::string("data.gel");
    config["is_restart"] = false;
    config["global_seed"] = 0;
    config["current_step"] = 0;
	config["log"] = json({
		{"name", "log.activerun" },
		{"flush", 100000 }
	});

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
		{"runtime_restart", true},
		{"runtime_restart_step", 10000000},
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
        if (restart.input_name.substr(restart.input_name.size() - 4, 4) == "json") {
            config.read_file(restart.input_name.c_str());
            config["input_file"] = std::string(restart.input_name);
        }
        else {
            return read_legacy_input(restart.input_name.c_str(), config);
        }

        config["is_restart"] = true;
        config["data_file"] = restart.data_name;
        config["current_step"] = restart.current_step;

        config["dump"]["dump_file"] = restart.output_name;
        if (!restart.thermo_name.empty()) {
            config["thermo"]["using_thermo"] = true;
            config["thermo"]["thermo_file"] = restart.thermo_name;
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
#ifdef READ_LEGACY_INPUT
        if (argc >= 3) {
            config["data_file"] = std::string(argv[2]);
        }
        return read_legacy_input(argv[1], config);
#else
        return 1;
#endif 
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

	json log_param = config["log"];
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


	/* Logger */
	MultiDumper logger_body({ stdout });
	if (log_param["name"] != "") {
		std::string log_name = log_param["name"];
		logger_body.add_file(log_name.c_str(), "w");
	}
	int flush_step = log_param["flush"];
	logger = &logger_body;


    /* Global configure variables */

    bool is_restart = config["is_restart"];
    std::string data_file = config["data_file"];
	logger->write_all(is_restart ? "Restarting...\n" : "Start new running task\n");

    /* Datafile input */

    DataFile datafile;
    datafile.read_data(data_file.c_str());

    System system;
    system.init(datafile);

    State state_init;
    state_init.init(datafile);

    Context context;

    /* group swim type 1 */

    bool has_swim = std::count(system.atom_type.begin(), system.atom_type.end(), 1) > 0;
    std::vector<int> group_swim(system.atom_num, 0);
    std::vector<int> group_passive(system.atom_num, 0);
    for (size_t i = 0; i < system.atom_num; i++) {
        group_swim[i] = system.atom_type[i] == 1 ? 1 : 0;
        group_passive[i] = group_swim[i] ? 0 : 1;
    }

    /* Forces */

    context.init_force_buffer(system, 3); // 3 forces

    /* fix 0 all langevin [BronwianForce.temp] [BrownianForce.temp] 0.0 */

    BrownianForce force_brown;
	force_brown.init(config.get_dict("BrownianForce"), system);
	force_brown.group_cache = &group_passive[0];
    
    /* fix 1 swim active [SwimForce.temp] [SwimForce.temp] 0.0 */

#ifndef THREE_DIMENSION
    SwimForce force_swim;
#else
    SwimForce3d force_swim;
#endif // !THREE_DIMENSION
    size_t swim_start = swim_param["swim_start"];
	if (has_swim) {
        force_swim.init(fix_swim_param, system);
        force_swim.group_cache = &group_swim[0];
	}

    double max_fix_cutoff = vec_max(system.get_attr("size"));

    /* */

    MorseForce force_morse;
	force_morse.init(pair_param, system);
    force_morse.group_cache = &group_swim[0];
    double max_pair_cutoff = force_morse.max_cutoff();

    /* Acceleration */

    double cell_size = util_param["cell_size"];
    int np = util_param["np"];
    int neighlist_step = util_param["neighlist_step"];
	size_t step_begin = config["current_step"];

	context.current_step = step_begin;
	context.init_timestep(time_param);
    context.init_pbc(system, true);
    context.init_neighlist(system.box, cell_size * std::max(max_pair_cutoff, max_fix_cutoff));
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

    /* dump 1 all custom [dump_step] [dump_file] id mol type x y z */

    std::string dump_file = dump_param["dump_file"].get<std::string>();
    size_t dump_step = dump_param["dump_step"].get<size_t>();

    TrajDumper trajdumper(dump_file.c_str(), is_restart);

    /* compute p_swim swim pressure 1 */

    FixPressureComputer p_swim;
    if (has_swim) {
        p_swim.init();
        p_swim.group_cache = &group_swim[0];
    }
#ifndef PRESSURE_BREAKDOWN
    /* compute p_morse all pressure pair */

    PressureComputer p_morse;
    p_morse.init(true);
#else
    std::vector<double> morse_pressure(3, 0.0), morse_energy(3, 0.0);
#endif

    /* thermo [thermo.step] */

    bool using_thermo = thermo_param["using_thermo"];
    std::string thermo_file = thermo_param["thermo_file"].get<std::string>();   // thermo output file
    bool do_thermo_sample = false;
    bool do_thermo_output = false;

    ThermoCounter thermocounter;
    thermocounter.init(thermo_param);
	if (step_begin > thermocounter.start_step()) {
		try {
			p_swim.read_restart("restart.pressure", context);
			context.pbc.update_image(state_init.pos);
			logger->write_all("Reading pressure restart succeeded\n");
		}
		catch (const std::runtime_error&) {
			thermocounter.check_start(step_begin);			
		}
	}

    std::vector<double> thermo_buffer;

#ifdef PRESSURE_BREAKDOWN
    ThermoDumper thermodumper(thermo_file.c_str(), *logger, { "P_Kinetics", "P_Swim", "P_Morse_pp", "P_Morse_aa", "P_Morse", "PE_pp", "PE_aa", "PE", "kT_Brown" }, is_restart);
    thermo_buffer.resize(9);
#else
    /* thermo_style custom step p_brown p_swim p_morse temp */

    LineDumper thermodumper(thermo_file.c_str(), { "P_Kinetics", "P_Swim", "P_Morse", "PE", "kT_Brown" }, true, is_restart);
    thermo_buffer.resize(5);
#endif

	/* Restart */
	RestartFile restart;

	bool write_restart = restart_param["write_restart"];
	bool write_runtime_restart = restart_param["runtime_restart"];
	int restart_step = restart_param["runtime_restart_step"];
	std::string input_file = config["input_file"];
	std::string restart_file = restart_param["restart_file"];
	std::string restart_data = restart_param["restart_data"];
	std::string restart_pressure = "restart.pressure";

	if (write_restart) {

		restart.input_name = input_file;
		restart.data_name = restart_data;
		restart.output_name = dump_file;
		restart.thermo_name = thermo_file;
		restart.current_step = step_begin;
	}

    /* timer normal every [sample_step] */

    Timer timer;
    timer.init({ "NeighList", "Force", "Integrator" });
    int timer_step = timer_param["sample_step"];



    // random seed

    int global_seed = config["global_seed"];
    if (global_seed == -1)global_seed = time(0);
	logger->write_all("Global seed=%d\n", global_seed);
	set_random_seed(global_seed);

    // building force cache

    integrator.update_cache(system, context);
    force_brown.update_cache(system, context);
    if (has_swim) {
		force_swim.update_cache(system, context);
    }
    force_morse.update_cache(system, context);
	p_swim.update_cache(system);

    /* run [time.step] */

	logger->write_all("\nRun %zd steps with timestep of %.4g...\n\n", context.total_steps - step_begin, context.timestep);

    State state = state_init;

    if (!is_restart) {

        trajdumper.dump(system, state, 0);

        if (using_thermo) {
	        thermodumper.dump_head();
        }
    }

	logger->flush_all();

    for (context.current_step = step_begin; context.current_step < context.total_steps; context.current_step++) {

        if (context.current_step % timer_step == 0) timer.set_checkpoint(-1);

        // update neighbourlist and pbc

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

        // test if it's a thermo step (then need to calculate energy, pe)

        if (using_thermo) {
            if (thermocounter.is_thermo_init_step(context.current_step)) {
                context.pbc.reset_location();
                context.pbc.update_image(state.pos);

				if (has_swim)p_swim.init_computing(context);
#ifndef PRESSURE_BREAKDOWN
                p_morse.update_cache(system, context);
#endif
            }

			// prevent first thermo output (which might be unstable)
			if (context.current_step - step_begin > thermocounter.step() / 2) {
				do_thermo_sample = thermocounter.is_thermo_sample_step(context.current_step);
				do_thermo_output = thermocounter.is_thermo_output_step(context.current_step);
			}
        }
        
        // update force
        
        if (context.current_step % timer_step == 0) timer.set_checkpoint(0);
        context.clear_buffer();

        force_morse.update_ahead(do_thermo_sample);
#ifdef DEBUG_FORCE
        try {
            force_morse.mp_update(context.pool, state, context);
        }
        catch (const std::runtime_error&) {
            trajdumper.dump(system, state, context.current_step / dump_step * dump_step);
            dump_snapshot(state, context);
            return 1;
        }
#else
		force_morse.mp_update(context.pool, state, context);    // morse force must update first
#endif
        force_brown.mp_update(context.pool, state, context.force_buffer[0]);
        if (has_swim && context.current_step >= swim_start) {
			force_swim.mp_update(context.pool, state, context.force_buffer[1]);
        }
		context.pool.wait();
		force_morse.update_later(context.force_buffer[2]);

        // compute

        if (do_thermo_sample) {
            context.pbc.update_image(state.pos);
            double temperature = force_brown.compute_temperature(context.force_buffer[0]);

            thermo_buffer[0] += state.pos.size() * temperature / system.volume();
            if (has_swim)thermo_buffer[1] += p_swim.compute_pressure(context, 1);
#ifndef PRESSURE_BREAKDOWN
            thermo_buffer[2] += p_morse.compute_pressure(context, state, 2);
            thermo_buffer[3] += force_morse.potential_energy();
            thermo_buffer[4] += temperature;
#else
            vec_div(force_morse.pressure, morse_pressure, (double) (DIMENSION * system.volume()));
            morse_energy = force_morse.energy;

            thermo_buffer[2] += morse_pressure[0];
            thermo_buffer[3] += morse_pressure[1];
            thermo_buffer[4] += vec_sum(morse_pressure);
            thermo_buffer[5] += morse_energy[0];
            thermo_buffer[6] += morse_energy[1];
            thermo_buffer[7] += vec_sum(morse_energy);
            thermo_buffer[8] += temperature;
#endif // !PRESSURE_BREAKDOWN
        }  

        // integrate

        if (context.current_step % timer_step == 0) timer.set_checkpoint(1);
        integrator.update(state, context);

        if (context.current_step % timer_step == 0) timer.set_checkpoint(2);

        // dump

        if ((context.current_step + 1) % dump_step == 0) {
            trajdumper.dump(system, state, context.current_step + 1);
        }

        // thermo

        if (do_thermo_output) {
            vec_div(thermo_buffer, thermo_buffer, (double)thermocounter.sample_count());
            thermodumper.dump(thermo_buffer, thermocounter.last_thermo_step(context.current_step + 1));
            vec_reset(thermo_buffer);
        }
           
		// flush
		if ((context.current_step + 1) % flush_step == 0) {
			logger->flush_all();
			trajdumper.flush();
		}

		if ((context.current_step + 1) % restart_step == 0 && write_runtime_restart) {

			state.write_data(datafile);
			datafile.write_data(restart.data_name.c_str());

			restart.current_step = context.current_step + 1;
			restart.write_restart(restart_file.c_str());
			p_swim.write_restart(restart_pressure.c_str(), context);

			logger->write_all("Restart written to %s\n", restart_file.c_str());
		}

    }

    timer.print((context.current_step - step_begin) / timer_step, context.current_step - step_begin);

    /* write_restart [restart_file] */
	if (write_restart) {
		state.write_data(datafile);
		datafile.write_data(restart.data_name.c_str());

		restart.current_step = context.current_step;
		restart.write_restart(restart_file.c_str());
		p_swim.write_restart(restart_pressure.c_str(), context);
		logger->write_all("\nRestart written to %s\n", restart_file.c_str());
	}


	return 0;
}
