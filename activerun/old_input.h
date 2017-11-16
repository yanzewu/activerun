#pragma once

/* Handling Legacy Input File (*.gel) */

#ifndef  ACTIVERUN_OLD_INPUT_H
#define ACTIVERUN_OLD_INPUT_H


#include "includes.h"
#include "Serializer.h"

struct InputParameter {
    double kT; //thermal energy
    double dt; //timestep normalized by shortest relaxation time
    double stepstR; //simulation duration in terms of tauR
    double zeta; //translational drag factor
    double overlap; //maximum distance particles can overlap (units of a
    double PeR; //rotational Peclet number = a/(UtauR
    int brownrot; //set to integer != 0 to include Brownian Rotation
    int browntrans;
    double PeS; //specify if brownrot == 0
    int attraction; //set to integer != 0 to include AO or Morse attraction
    int swimattract; //set to integer != 0 for swimmers to also feel attraction
    int morse;//set to integer != 0 to use the Morse potential for attraction+repulsion/repulsion (if attraction == 0
    int HSrepulsion; //set to integer != 0 to include HS or Morse repulsion
    double Rg; //sets the cutoff of AO or Morse attraction
    double kappa; //Morse screening length (units of 1/a
    double Um; //sets the well depth of the AO or Morse attraction
    int output; //output frequency (units of timestep
    int swimstart; //timestep that active colloids begin to swim
    int spstart; //timestep that swim pressure measurement begins
    double mincellL;
    double U;
    double viscosity;
    int np;

    int read_input(const char* filename);
};

// read input file
int InputParameter::read_input(const char* filename) {
    FILE* file_params;
    if ((file_params = fopen("in.gel", "r")) == NULL) {
        printf("%s file could not be opened\n", filename);
        //throw;
        return 1;
    }

    fscanf(file_params, "%*s%lf", &kT); //thermal energy
    fscanf(file_params, "%*s%lf", &dt); //timestep normalized by shortest relaxation time
    fscanf(file_params, "%*s%lf", &stepstR); //simulation duration in terms of tauR
    fscanf(file_params, "%*s%lf", &zeta); //translational drag factor
    fscanf(file_params, "%*s%lf", &overlap); //maximum distance particles can overlap (units of a)
    fscanf(file_params, "%*s%lf", &PeR); //rotational Peclet number = a/(UtauR)
    fscanf(file_params, "%*s%d", &brownrot); //set to integer != 0 to include Brownian Rotation
    fscanf(file_params, "%*s%i", &browntrans); //set to integer != 0 to also include Brownian translation for the swimmers
    fscanf(file_params, "%*s%lf", &PeS); //specify if brownrot == 0
    fscanf(file_params, "%*s%d", &attraction); //set to integer != 0 to include AO or Morse attraction
    fscanf(file_params, "%*s%d", &swimattract); //set to integer != 0 for swimmers to also feel attraction
    fscanf(file_params, "%*s%d", &morse);//set to integer != 0 to use the Morse potential for attraction+repulsion/repulsion (if attraction == 0)
    fscanf(file_params, "%*s%i", &HSrepulsion);//set to integer !=0 to use HS repulsion (default with AO potential) in conjunction with Morse Attraction
    fscanf(file_params, "%*s%lf", &Rg); //sets the cutoff of AO or Morse attraction
    fscanf(file_params, "%*s%lf", &kappa); //Morse screening length (units of 1/a)
    fscanf(file_params, "%*s%lf", &Um); //sets the well depth of the AO or Morse attraction
    fscanf(file_params, "%*s%d", &output); //output frequency (units of timestep)
    fscanf(file_params, "%*s%d", &swimstart); //timestep that active colloids begin to swim
    fscanf(file_params, "%*s%d", &spstart); //timestep that swim pressure measurement begins
    fscanf(file_params, "%*s%lf", &mincellL); //length of unit cell in cell list normalized by largest interaction length
    fclose(file_params);

    np = 4;

    if (attraction != 0 && morse == 0)
        HSrepulsion = 1;

    viscosity = zeta / 6 / M_PI / 2.0;

    return 0;
}

int read_legacy_input(const char* filename, Serializer& config) {

    InputParameter input;
    if (input.read_input(filename)) {
        return 1;
    }

    config["input_file"] = std::string(filename);
    config["global_seed"] = -1;

    config["BrownianForce"]["temp"] = input.kT;
    config["BrownianForce"]["neta"] = input.viscosity;

    config["SwimForce"]["swim_temp"] = input.kT;
    config["SwimForce"]["neta"] = input.viscosity;
    config["SwimForce"]["brownian"] = (input.brownrot == 1.0);
    config["SwimForce"]["PeR"] = input.PeR;
    config["SwimForce"]["swim_start"] = input.spstart;

    config["MorseForce"]["kappa"] = input.kappa;
    config["MorseForce"]["Um"] = input.Um;
    config["MorseForce"]["cutoff"] = input.Rg * 2 + 1.0;

    config["run"]["timestep"] = input.dt;
    config["run"]["steps"] = input.stepstR / input.dt;

    config["dump"]["dump_step"] = input.output;
    config["thermo"]["thermo_step"] = input.output;
    config["thermo"]["thermo_start"] = input.spstart;

    config["util"]["cell_size"] = input.mincellL;

    return 0;
}

#endif // ! ACTIVERUN_OLD_INPUT_H
