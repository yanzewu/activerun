#pragma once

#ifndef ACTIVERUN_INPUT_H
#define ACTIVERUN_INPUT_H


#include "includes.h"

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

// LAMMPS Data File
struct DataFile {
public:
	int atom_num;
	int atom_type_num;

	std::vector<double> type_mass;
	std::vector<std::vector<double> > type_pair_coeff;

	std::vector<Vec3> pos;
	std::vector<int> atom_type;
	std::vector<int> atom_group;

	Vec3 box, box_start;

	void init_size();

	// Read LAMMPS Data File
	int read_data(const char*);

	int write_data(const char*);
};


#endif // !ACTIVERUN_INPUT_H

