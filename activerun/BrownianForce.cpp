#include "Force.h"


BrownianForce::BrownianForce()
{
	this->is_direct = false;
	this->is_paired = false;
	this->is_potential_force = false;
}

BrownianForce::~BrownianForce() {
}

void BrownianForce::init(const Dict& params, System& system) {

	logger->write_all("\nInitializing bronwian force\n\n");

	temperature = params.get("temp", 1.0);
	logger->write_all("Temperature=%.5g\n", temperature);

	force_coeff_cache.resize(system.atom_num);
	random_cache.resize(system.atom_num);

	// zeta = 6 pi eta r^2
	if (system.attribute_names.find("zeta") == system.attribute_names.end()) {
		logger->write_all("Damp cache not found, building...\n");
		logger->write_all("Viscosity=%.5g\n", params.get("neta", 1.0));

		auto& zeta = system.add_attr("zeta");
		system.set_name("size", 2);

		for (size_t i = 0; i < system.atom_num; i++) {
			zeta[i] = 6 * M_PI * params.get("neta", 1.0) * system.atom_attributes[2][i];
		}
	}
	else {
		logger->write_all("Using cached damp\n");
	}
}

void BrownianForce::init_mpi(int thread_count) {
	logger->write_all("Brownian force: ");
	using_thread = thread_count > 0;
	logger->write_all(using_thread ? "Using thread pool\n" : "Using main thread\n");
}

void BrownianForce::update_ahead(State& state, std::vector<Vec>& force_buffer) {

}

// update
void BrownianForce::update(const State& state, std::vector<Vec>& force_buffer) {

	
	for (auto& rc : random_cache) {	/* If you want to be consistent with history, reverse the order */
#ifdef THREE_DIMENSION
		rc[2] = rand_uniform() - 0.5;
#endif // THREE_DIMENSION
		rc[1] = rand_uniform() - 0.5;
		rc[0] = rand_uniform() - 0.5;
	}

	for (size_t i = 0; i < force_buffer.size(); ++i) {
		if (!group_cache[i]) continue;
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

double BrownianForce::compute_temperature(const std::vector<Vec>& force)const
{
    double T = 0.0;
    for (size_t i = 0; i < force.size(); i++) {
		if (!group_cache[i]) continue;
        T += (force[i] / force_coeff_cache[i]).norm2();
    }
	return T * 12 / DIMENSION / std::count(group_cache, group_cache + force.size(), 1) * temperature;
}

