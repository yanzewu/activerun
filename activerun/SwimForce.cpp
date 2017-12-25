#include "Force.h"

SwimForce::SwimForce() 
{
	this->is_direct = false;
	this->is_paired = false;
	this->is_potential_force = false;
}

void SwimForce::init(const Dict& params, System& system) {

	printf("\nInitializing Swim Force\n\n");

	try {
		Pe_R.resize(system.atom_num, params.at("PeR"));
	}
	catch (const std::out_of_range&) {
		printf("Error: PeR not found\n");
		throw;
	}
    printf("Rotation Peclet=%.4f\n", Pe_R[0]);

	temperature = params.get("swim_temp", 1.0);
	brownian_rotation = (bool)params.get("brownian", 1.0);
	temperature = params.get("swim_temp", 1.0);
	if (brownian_rotation) {
		printf("Using brownian rotation\nrotation temperature=%.4f\n", temperature);
	}
	else {
		Pe_S.resize(system.atom_num, params.get("PeS", 1.0));
		printf("Using self-defined rotation\nSwim Peclet=%.4f\n", params.get("PeS", 1.0));
	}

	my_type = (int)params.get("type", 1.0);
	printf("Swim atom type: %d\n", my_type);

	int rand_seed = (int)params.get("init_seed", 0.0);
	printf("Randomize initial configuration...\nseed=%d\n", rand_seed);
	srand(rand_seed);

	angle_cache.resize(system.atom_num);
	for (auto& angle : angle_cache) {
		angle = rand_uniform() * 2 * M_PI;
	}

	torque_coeff_cache.resize(system.atom_num);
	force_coeff_cache.resize(system.atom_num);
	angular_momentum_cache.resize(system.atom_num);
	rot_coeff_cache.resize(system.atom_num);

	// zeta = 6 pi eta r^2
	if (system.attribute_names.find("zeta") == system.attribute_names.end()) {
		printf("Damp cache not found, building...\n");
		printf("Viscosity=%.5g\n", params.get("neta", 1.0));

		auto& zeta = system.add_attr("zeta");
		system.set_name("size", 2);

		for (size_t i = 0; i < system.atom_num; i++) {
			zeta[i] = 6 * M_PI * params.get("neta", 1.0) * system.atom_attributes[2][i];
		}
	}
	else {
		printf("Using cached damp\n");
	}
}

void SwimForce::init_mpi(int thread_count) {
	printf("Swim force: ");
	using_thread = thread_count > 0;
	printf(using_thread ? "Using thread pool\n" : "Using main thread\n");
}

void SwimForce::update_ahead(State& state, std::vector<Vec>& F) {

}

void SwimForce::update(const State& state, std::vector<Vec>& force_buffer) {
	for (int i = 0; i < force_buffer.size(); ++i) {
		if (!group_cache[i]) continue;
		angular_momentum_cache[i] = torque_coeff_cache[i] * (rand_uniform() - 0.5);
	}

	for (int i = 0; i < force_buffer.size(); i++) {
		angle_cache[i] += angular_momentum_cache[i] * rot_coeff_cache[i];
	}

	for (int i = 0; i < force_buffer.size(); i++) {
		if (!group_cache[i])continue;
#ifdef THREE_DIMENSION
		force_buffer[i] = Vec3(cos(angle_cache[i]), sin(angle_cache[i]), 0) * force_coeff_cache[i];
#else
		force_buffer[i] = Vec2(cos(angle_cache[i]), sin(angle_cache[i])) * force_coeff_cache[i];

#endif // THREE_DIMENSION
	}

}

void SwimForce::mp_update(FixedThreadPool& pool, const State& state, std::vector<Vec>& force_buffer) {
	if (using_thread) {
		pool.submit(std::bind(&SwimForce::update, this, std::ref(state), std::ref(force_buffer)));
	}
	else {
		update(state, force_buffer);
	}
}

void SwimForce::update_cache(const System& system, const Context& context) {

	const std::vector<double>& atom_size = system.get_attr("size");
	const std::vector<double>& zeta = system.get_attr("zeta");
	std::vector<double> zeta_R(system.atom_num, 0.0);
	std::vector<double> tau_R(system.atom_num, 0.0);

	// zetaR = zeta*D^2/3
	for (size_t i = 0; i < zeta.size(); i++) {
		zeta_R[i] = zeta[i] * atom_size[i] * atom_size[i] / 3;
	}

	if (!brownian_rotation) {
		// tauR = 0.75 zetaR/(T*PeR*PeS)
		for (size_t i = 0; i < tau_R.size(); i++) {
			tau_R[i] = zeta_R[i] / (temperature * Pe_R[i] * Pe_S[i]) * 0.75;
		}
	}
	else {
		// tauR = zetaR/kT
		for (size_t i = 0; i < tau_R.size(); i++) {
			tau_R[i] = zeta_R[i] / temperature;
		}
	}

	for (size_t i = 0; i < zeta_R.size(); i++) {
		rot_coeff_cache[i] = context.timestep / zeta_R[i];
	}

	// swim speed
	for (int i = 0; i < system.atom_num; i++) {
		force_coeff_cache[i] = 0.5 * atom_size[i] / (Pe_R[i] * tau_R[i]) * zeta[i];
	}

	for (int i = 0; i < system.atom_num; i++) {
		torque_coeff_cache[i] = sqrt(
			24.0 * zeta_R[i] * zeta_R[i] / tau_R[i] / context.timestep);
	}
}
