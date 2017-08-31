# ActiveRun Usage

## Usage

    activerun [legacy input] [legacy data]
    activerun -i [json input]
    activerun -r [restart file]

- Legacy Input: Row-specified Data
- Json Input: Input with json format. A json with default value will be created if input not found.


## Sample Input File

### Legacy Input

    kT 1.0
    dt 0.000005             # absolute timestep
    stepstauB 1             # step for tauB (not tauR)
    zeta(stokes) 1          # zeta when radius=1.0
    overlap(a) 0.0000001
    PeR 0.05
    brownrot 1              # not work
    browntrans 1            # not work
    PeS 0.0
    attraction 1            # not work
    swimattract 1           # not work
    morse 1
    HS 0                    # not work
    Rg(a) 0.1
    kappa 30
    Um(kT) 5 
    output 10000
    swimstart 0
    spstart 0
    mincellL 2

Currently the type of forces need to be changed manually.

### Json Input

    {
        "BrownianForce":{
            "neta":0.026525823848649224,    # viscosity coeff
            "temp":1.0                      # temperature for Brownian Force
        },
        "MorseForce":{
            "Um":5.0,
            "cutoff":1.20,
            "kappa":30.0
        },
        "SwimForce":{
            "PeR":0.01,
            "brownian":true,                # not work
            "init_seed":0,                  # seed for randomize angle
            "neta":0.026525823848649224,
            "swim_start":0,
            "swim_temp":1.0
        },
        "current_step":0,
        "data_file":"data.gel",
        "dump":{
            "dump_file":"outputconfig.lammpstrj",
            "dump_step":1000
        },
        "global_seed":0,                    # random seed for whole running session
        "is_restart":false,
        "restart":{                         # configure for restart
            "restart_data":"restart.data",
            "restart_file":"restart.rst",
            "write_restart":true
        },
        "run":{
            "step":200000,                  # if start from restart file, it's the total step
            "timestep":5e-06                # absolute timestep
        },
        "thermo":{                          # thermodynamic property (currently only pressure and temperature)
            "average_range":10,             # average so many step around sample point
            "average_step":1,
            "compute_temp":true,
            "thermo_file":"pressure.output",
            "thermo_start":0,
            "thermo_step":100,              # calculate step
            "using_thermo":true
        },
        "util":{
            "cell_size":2.0,
            "neighlist_step":1,
            "np":2                          # core number 
        }
    }
