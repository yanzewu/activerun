
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

    /* Global configure variables */

    bool is_restart = config["is_restart"];
    std::string data_file = config["data_file"];
    size_t step_begin = config["current_step"];

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
    double max_pair_cutoff = force_morse.max_cutoff();

    /* Acceleration */

    double cell_size = util_param["cell_size"];
    int np = util_param["np"];
    int neighlist_step = util_param["neighlist_step"];

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

    /* compute p_brown all pressure 0 */
    /* compute p_swim swim pressure 1 */

    FixPressureComputer p_brown, p_swim;
    p_brown.init();
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
    
    ThermoCounter thermocounter;
    thermocounter.init(thermo_param);
    thermocounter.check_start(step_begin);

#ifdef PRESSURE_BREAKDOWN
    LineDumper thermodumper(thermo_file.c_str(), { "P_Brown", "P_Swim", "P_Morse", "P_Morse_pp", "P_Morse_aa", "P_Morse", "PE", "PE_pp", "PE_aa", "T" }, true, is_restart);
#else
    /* thermo_style custom step p_brown p_swim p_morse temp */

    LineDumper thermodumper(thermo_file.c_str(), { "P_Brown", "P_Swim", "P_Morse", "PE", "T" }, true, is_restart);
    std::vector<double> thermo_buffer(5, 0.0);
#endif

    /* timer normal every [sample_step] */

    Timer timer;
    timer.init({ "NeighList", "Force", "Integrator" });
    int timer_step = timer_param["sample_step"];

    bool do_thermo_sample = false;
    bool do_thermo_output = false;

    // random seed

    int global_seed = config["global_seed"];
	set_random_seed(global_seed);

    // building force cache

    integrator.update_cache(system, context);
    force_brown.update_cache(system, context);
    if (has_swim) {
		force_swim.update_cache(system, context);
    }
    force_morse.update_cache(system, context);

    /* run [time.step] */

	printf("\nStart running session...\n\n");

    State state = state_init;

    if (!is_restart) {

        trajdumper.dump(system, state, 0);

        if (using_thermo) {
	        thermodumper.dump_head();
        }
    }

    for (context.current_step = step_begin; context.current_step < context.total_steps; context.current_step++) {

        // test if it's a thermo step (then need to calculate energy, pe)

        if (using_thermo) {
            if (thermocounter.is_thermo_init_step(context.current_step)) {
                context.pbc.reset_location();
                context.pbc.update_image(state.pos);

                p_brown.update_cache(system, context);
                if (has_swim)p_swim.update_cache(system, context);
#ifndef PRESSURE_BREAKDOWN
                p_morse.update_cache(system, context);
#endif
            }
            do_thermo_sample = thermocounter.is_thermo_sample_step(context.current_step + 1);
            do_thermo_output = thermocounter.is_thermo_output_step(context.current_step + 1);
        }

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

        // update force
        
        if (context.current_step % timer_step == 0) timer.set_checkpoint(0);
        context.clear_buffer();

        //try {
        force_morse.update_ahead(do_thermo_sample);
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

        // integrate

        if (context.current_step % timer_step == 0) timer.set_checkpoint(1);
        integrator.update(state, context);

        if (context.current_step % timer_step == 0) timer.set_checkpoint(2);

        // dump

        if ((context.current_step + 1) % dump_step == 0) {
            trajdumper.dump(system, state, context.current_step + 1);
        }

        // compute

        if (do_thermo_sample) {
            context.pbc.update_image(state.pos);
            thermo_buffer[0] += p_brown.compute_pressure(context, 0);
            if (has_swim)thermo_buffer[1] += p_swim.compute_pressure(context, 1);
#ifndef PRESSURE_BREAKDOWN
            thermo_buffer[2] += p_morse.compute_pressure(context, state, 2);
            thermo_buffer[3] += force_morse.potential_energy();
#else
            for (int i = 0; i < 3; i++)morse_pressure[i] += force_morse.pressure[i] / (DIMENSION * system.volume());
            for (int i = 0; i < 3; i++)morse_energy[i] += force_morse.energy[i];
#endif // !PRESSURE_BREAKDOWN
            thermo_buffer[4] += integrator.update_temperature();
        }

        // thermo

        if (do_thermo_output) {
#ifndef PRESSURE_BREAKDOWN
            vec_div(thermo_buffer, thermo_buffer, (double)thermocounter.sample_count());
            thermodumper.dump(thermo_buffer, thermocounter.last_thermo_step(context.current_step + 1));
#else
            thermodumper.dump({
                pressure[0] / thermocounter.sample_count(),
                pressure[1] / thermocounter.sample_count(),
                morse_pressure[0] / thermocounter.sample_count(),
                morse_pressure[1] / thermocounter.sample_count(),
                (morse_pressure[0] + morse_pressure[1] + morse_pressure[2]) / thermocounter.sample_count(),
                morse_energy[2] / thermocounter.sample_count(),
                morse_energy[2] / thermocounter.sample_count(),
                (morse_energy[0] + morse_energy[1] + morse_energy[2]) / thermocounter.sample_count(),
                temperature / thermocounter.sample_count() },
                thermocounter.last_thermo_step(context.current_step + 1));
            memset(&morse_pressure[0], 0, 3 * sizeof(double));
            memset(&morse_energy[0], 0, 3 * sizeof(double));
#endif 
            vec_reset(thermo_buffer);
        }
            
    }

    timer.print((context.current_step - step_begin) / timer_step, context.current_step - step_begin);

    /* write_restart [restart_file] */

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