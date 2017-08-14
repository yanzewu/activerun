#include "Input.h"


char* rstrip(char* str, char c) {
	while (str[strlen(str) - 1] == c) {
		str[strlen(str) - 1] = '\0';
	}
	return str;
}


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

	if (attraction != 0 && morse == 0)
		HSrepulsion = 1;

	viscosity = zeta / 6 / M_PI / 2.0;

	return 0;
}


void DataFile::init_size() {
	type_mass.resize(atom_type_num, 0.0);
	type_pair_coeff.resize(atom_type_num, { 0.0, 0.0, 0.0, 0.0 });

	pos.resize(atom_num);
	atom_type.resize(atom_num);
	atom_group.resize(atom_num);
}

/* Read LAMMPS Data File */
int DataFile::read_data(const char* filename) {

	FILE *file_data;
	if ((file_data = fopen(filename, "r")) == NULL) {
		fprintf(stderr, "%s could not be opened\n", filename);
		return 1;
	}

	const int buffer_size = 128;
	char buffer[buffer_size];

	// head section
	fgets(buffer, buffer_size, file_data);
	if (strcmp(buffer, "LAMMPS Data File\n")) {
		fprintf(stderr, "Unrecognized file_data file!\n");
		return 1;
	}
	fgets(buffer, buffer_size, file_data);
	char name_buffer[128];

	// num section

	while (true) {
		fgets(buffer, buffer_size, file_data);
		if (buffer[0] == '\n')break; // empty

		int tmp;
		sscanf(buffer, "%d%s\n", &tmp, name_buffer);

		if (strcmp(name_buffer, "atoms") == 0) {
			atom_num = tmp;
		}
		else {

		}
	}

	// type num section

	while (true) {
		fgets(buffer, buffer_size, file_data);
		if (buffer[0] == '\n')break; // empty

		int tmp;
		sscanf(buffer, "%d%s", &tmp, name_buffer);

		if (strcmp(name_buffer, "atom") == 0) {
			atom_type_num = tmp;
		}
		else {

		}
	}

	init_size();

	// box section

	for (int i = 0; i < 2; i++) {
		fgets(buffer, buffer_size, file_data);
		if (buffer[0] == '\n') {
			fprintf(stderr, "Invalid box size section\n");
			return 1;
		};

		sscanf(buffer, "%lf%lf", &box_start[i], &box[i]);
	}
	while (true) {
		fgets(buffer, buffer_size, file_data);
		if (buffer[0] == '\n')break; // empty
	}


	// type coeff section

	enum SectionType {
		NONE, MASS, PAIR_COEFF, ATOMS
	} section_type;

	while (true) {
		if (!fgets(buffer, buffer_size, file_data)) {
			break;
		}

		rstrip(buffer, '\n');

		if (strcmp(buffer, "Masses") == 0) {
			section_type = SectionType::MASS;
		}
		else if (strcmp(buffer, "Pair Coeffs") == 0) {
			section_type = SectionType::PAIR_COEFF;
		}
		else if (strcmp(buffer, "Atoms") == 0) {
			section_type = SectionType::ATOMS;
		}
		else {
			section_type = SectionType::NONE;
		}

		strcpy(name_buffer, buffer);
		fgets(buffer, buffer_size, file_data); // empty line

		switch (section_type)
		{
		case NONE:
			fprintf(stderr, "Unknown section %s\n", name_buffer);
			return 1;
		case MASS:
			for (int i = 0; i < atom_type_num; i++) {
				int type;
				fgets(buffer, buffer_size, file_data);
				sscanf(buffer, "%d", &type);
				sscanf(buffer, "%*s%lf", &type_mass[type - 1]);
			}
			break;
		case PAIR_COEFF:
			for (int i = 0; i < atom_type_num; i++) {
				int type;
				fgets(buffer, buffer_size, file_data);
				sscanf(buffer, "%d", &type);
				sscanf(buffer, "%*s%lf%lf%lf%lf", &type_pair_coeff[type - 1][0],
					&type_pair_coeff[type - 1][1],
					&type_pair_coeff[type - 1][2],
					&type_pair_coeff[type - 1][3]);
			}
			break;
		case ATOMS:
			for (int j = 0; j < atom_num; ++j) {

				int id;

				fgets(buffer, buffer_size, file_data);
				sscanf(buffer, "%d", &id);
				sscanf(buffer, "%*d%d%d%lf%lf%*i", &atom_group[id - 1],
					&atom_type[id - 1],
					&pos[id - 1][0],
					&pos[id - 1][1]);

			}
			break;
		default:
			break;
		}

		if (!fgets(buffer, buffer_size, file_data)) {
			break;
		}
	}
	fclose(file_data);
	return 0;
}
