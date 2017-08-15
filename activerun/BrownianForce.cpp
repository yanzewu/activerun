#include "Force.h"


BrownianForce::BrownianForce() : Force(false, false)
{

}

BrownianForce::~BrownianForce() {
}

void BrownianForce::init(const Dict& params, System& system) {

	printf("\nInitializing bronwian force\n\n");

	temperature = params.get("temp", 1.0);
	printf("Temperature=%.5g\n", temperature);

	group_cache.resize(system.atom_num);
	force_coeff_cache.resize(system.atom_num);
	random_cache.resize(system.atom_num);

	// zeta = 6 pi eta r^2
	if (system.attribute_names.find("zeta") == system.attribute_names.end()) {
		printf("Damp cache not found, buidling...\n");
		printf("Viscosity=%.5g\n", params.get("neta", 1.0));

		auto& zeta = system.add_attr("zeta");
		system.set_name("size", 2);

		for (size_t i = 0; i < system.atom_num; i++) {
			zeta[i] = 6 * M_PI * params.get("neta", 1.0) * system.atom_attributes[2][i];
		}
	}

}

void BrownianForce::init_mpi(int thread_count) {
	using_thread = thread_count > 0;
}

void BrownianForce::update_ahead(State& state, std::vector<Vec>& force_buffer) {

}

// update
void BrownianForce::update(const State& state, std::vector<Vec>& force_buffer) {

	for (auto& rc : random_cache) {
		rc[0] = rand_uniform() - 0.5;
		rc[1] = rand_uniform() - 0.5;
#ifdef THREE_DIMENSION
		rc[2] = rand_uniform() - 0.5;
#endif // THREE_DIMENSION

	}

	for (size_t i = 0; i < force_buffer.size(); ++i) {
		force_buffer[i] = random_cache[i] * force_coeff_cache[i];
	}
}

void BrownianForce::mp_update(FixedThreadPool& pool, const State& state, std::vector<Vec>& force_buffer) {
	if (using_thread) {
		pool.submit(std::bind(&BrownianForce::update, this, std::ref(state), std::ref(force_buffer)));
	}
	else {
		update(state, force_buffer);
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

