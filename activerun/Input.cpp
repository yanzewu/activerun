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

	np = 4;

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
		throw std::exception("IO Error");
	}

	const int buffer_size = 128;
	char buffer[buffer_size];

	// head section
	fgets(buffer, buffer_size, file_data);
	if (strcmp(buffer, "LAMMPS Data File\n")) {
		fprintf(stderr, "Unrecognized file_data file!\n");
		throw std::exception("IO Error");
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

int DataFile::write_data(const char * filename)
{
	FILE* file_data = fopen(filename, "w");
	fprintf(file_data, "LAMMPS Data File\n\n");

	fprintf(file_data, "%d atoms\n", atom_num);
	fprintf(file_data, "0 bonds\n0 angles\n0 dihedrals\n0 impropers\n");
	fprintf(file_data, "\n");

	fprintf(file_data, "%d atom types\n", atom_type_num);
	fprintf(file_data, "\n");

	fprintf(file_data, "%f %f xlo xhi\n", box_start[0], box[0]);
	fprintf(file_data, "%f %f ylo yhi\n", box_start[1], box[1]);
#ifdef THREE_DIMENSION
	fprintf(file_data, "%f %f zlo zhi\n", box_start[2], box[2]);
#else
	fprintf(file_data, "%f %f zlo zhi\n", -0.25, 0.25);
#endif // THREE_DIMENSION
	fprintf(file_data, "%f %f %f xy xz yz\n", 0.0, 0.0, 0.0);
	fprintf(file_data, "\n");

	fprintf(file_data, "Masses\n\n");
	for (int i = 0; i < atom_type_num; i++) {
		fprintf(file_data, "%d %f\n", i + 1, type_mass[i]);
	}
	fprintf(file_data, "\n");

	fprintf(file_data, "Pair Coeffs\n\n");
	for (int i = 0; i < atom_type_num; i++) {
		fprintf(file_data, "%d", i + 1);
		for (int j = 0; j < type_pair_coeff[i].size(); j++) {
			fprintf(file_data, " %f", type_pair_coeff[i][j]);
		}
		fprintf(file_data, "\n");
	}
	fprintf(file_data, "\n");

	fprintf(file_data, "Atoms\n\n");
	for (int i = 0; i < atom_num; i++) {
		fprintf(file_data, "%d %d %d", i + 1, atom_group[i], atom_type[i]);
#ifdef THREE_DIMENSION
		fprintf(file_data, " %f %f %f\n", pos[i][0], pos[i][1], pos[i][2]);
#else
		fprintf(file_data, " %f %f %f\n", pos[i][0], pos[i][1], 0);
#endif // THREE_DIMENSION
	}

	fclose(file_data);

	return 0;
}
