#include "Compute.h"
#include "output.h"

void FixPressureComputer::init() {
	group_cache = nullptr;
}

double FixPressureComputer::compute_pressure(const Context& context, int force_id) {
	double pressure = 0.0;

	for (size_t j = 0; j < context.force_buffer[force_id].size(); j++) {
		if (group_cache && !group_cache[j])continue;
		pressure += context.force_buffer[force_id][j].dot(context.pbc.image_cache[j] - init_pos[j]);
	}

	return pressure / (volume * DIMENSION);
}

void FixPressureComputer::update_cache(const System& system) {
	volume = system.volume();
}

void FixPressureComputer::init_computing(const Context& context) {
	init_pos = context.pbc.image_cache;
}

void FixPressureComputer::read_restart(const char* filename, Context& context) {
	FILE* ifile = fopen(filename, "r");
	init_pos.resize(context.pbc.location_cache.size());
	
	for (size_t i = 0; i < context.pbc.location_cache.size(); i++){
		for (int j = 0; j < DIMENSION; j++) {
			fscanf(ifile, "%le", &init_pos[i][j]);
		}
		for (int j = 0; j < DIMENSION; j++) {
			fscanf(ifile, "%d", &context.pbc.location_cache[i][j]);
		}
	}
	fclose(ifile);
}

void FixPressureComputer::write_restart(const char* filename, const Context& context)const {
	FILE* ofile = fopen(filename, "w");
	for (size_t i = 0; i < init_pos.size(); i++) {
		for (int j = 0; j < DIMENSION; j++) {
			fprintf(ofile, "%.18e ", init_pos[i][j]);
		}
		for (int j = 0; j < DIMENSION; j++) {
			fprintf(ofile, "%d ", context.pbc.location_cache[i][j]);
		}

		fprintf(ofile, "\n");
	}
	fclose(ofile);
}


void PressureComputer::init(double using_pbc) {
	this->using_pbc = using_pbc;
	group_cache = nullptr;
}

double PressureComputer::compute_pressure(const Context& context, const State& state, int force_id) {

	double pressure = 0.0;

	for (size_t j = 0; j < context.force_buffer[force_id].size(); j++) {
		if (group_cache && !group_cache[j])continue;
		pressure += context.force_buffer[force_id][j].dot(state.pos[j]);
	}

	std::vector<size_t> ghost;
	std::vector<Vecd> offset;

	context.neigh_list->compute_ghost_pos(ghost, offset);
	for (size_t i = 0; i < ghost.size(); i++) {
		if (group_cache && !group_cache[ghost[i]])continue;
		pressure += context.force_buffer[force_id][ghost[i]].dot(state.pos[ghost[i]] + offset[i].to_float() * box);
	}

	return pressure / (DIMENSION * volume);
}

void PressureComputer::update_cache(const System& system, const Context& context) {
	volume = system.volume();
	box = system.box;
}


void KineticEnergyComputer::init() {

}

double KineticEnergyComputer::compute_temperature(const std::vector<Vec>& velocity) {
	double temperature = 0.0;
	for (const auto& v : velocity) {
		temperature += v.norm2();
	}
	return temperature / (velocity.size() * DIMENSION);

}
