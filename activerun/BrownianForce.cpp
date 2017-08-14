#include "Force.h"


BrownianForce::BrownianForce() : Force(false, false)
{

}

BrownianForce::~BrownianForce() {
}

void BrownianForce::init(const InputParameter& input_params, System& system) {
	temperature = input_params.kT;

	group_cache.resize(system.atom_num);
	force_coeff_cache.resize(system.atom_num);
	random_cache.resize(system.atom_num);

	// zeta = 6 pi eta r^2
	if (system.attribute_names.find("zeta") == system.attribute_names.end()) {
		auto& zeta = system.add_attr("zeta");
		system.set_name("size", 2);

		for (size_t i = 0; i < system.atom_num; i++) {
			zeta[i] = 6 * M_PI * input_params.viscosity * system.atom_attributes[2][i];
		}
	}

}

void BrownianForce::update_ahead(State& state, std::vector<Vec2>& force_buffer) {

}

// update
void BrownianForce::update(const State& state, std::vector<Vec2>& force_buffer) {

	for (size_t i = 0; i < random_cache.size(); ++i) {
		random_cache[i] = Vec2(rand_uniform() - 0.5, rand_uniform() - 0.5);
	}

	for (size_t i = 0; i < force_buffer.size(); ++i) {
		force_buffer[i] = random_cache[i] * force_coeff_cache[i];
	}
}

void BrownianForce::update_cache(const System& system, const Context& context) {

	const auto& zeta = system.get_attr("zeta");
	for (size_t i = 0; i < system.atom_num; i++) {
		force_coeff_cache[i] = sqrt(24.0 * temperature * zeta[i] / context.timestep);
	}
}

double BrownianForce::compute_pressure(const State& state, const std::vector<Vec2>& force) {
	return 0.0;
}

double BrownianForce::compute_energy(const State& state) {
	return 2 * temperature;
}

