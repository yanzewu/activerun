#pragma once

#ifndef ACTIVERUN_INPUT_H
#define ACTIVERUN_INPUT_H


#include "includes.h"


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


struct RestartFile {
    char input_name[128];
    char data_name[128];
    char output_name[128];
    char thermo_name[128];
    size_t current_step;

    RestartFile() {
        thermo_name[0] = '\0';
    }

    int read_restart(const char*);

    void write_restart(const char*);
};

#endif // !ACTIVERUN_INPUT_H

