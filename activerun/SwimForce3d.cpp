#include "Force.h"

#ifdef THREE_DIMENSION

SwimForce3d::SwimForce3d()
{
	this->is_direct = false;
	this->is_paired = false;
	this->is_potential_force = false;
}

void SwimForce3d::init(const Dict& params, System& system) {

	printf("\nInitializing 3D Swim Force\n\n");

	try {
		Pe_R.resize(system.atom_num, params.at("PeR"));
	}
	catch (const std::out_of_range&) {
		printf("Error: PeR not found\n");
		throw;
	}
    printf("Rotation Peclet=%.4f", Pe_R[0]);

	brownian_rotation = (bool)params.get("brownian", 1.0);
	if (brownian_rotation) {
		temperature = params.get("swim_temp", 1.0);
		printf("Using brownian rotation\nrotation temperature=%.4f\n", temperature);
	}
	else {
		Pe_S.resize(system.atom_num, params.get("PeS", 1.0));
		printf("Using self-defined rotation\nSwim Peclet=%.4f\n", params.get("PeS", 1.0));
	}

	my_type = (int)params.get("type", 1.0);
	printf("Swim atom type: %d\n", my_type);
	group_cache.resize(system.atom_num);

	int rand_seed = (int)params.get("init_seed", 0.0);
	printf("Randomize initial configuration...\nseed=%d\n", rand_seed);
	srand(rand_seed);

	direct_cache.resize(system.atom_num);
	for (auto& direct : direct_cache) {
        direct[2] = rand_uniform() * 2 - 1.0;
        double theta = rand_uniform() * 2 * M_PI;
        direct[0] = cos(theta) * sqrt(1 - direct[2] * direct[2]);
        direct[1] = sin(theta) * sqrt(1 - direct[2] * direct[2]);
	}

	torque_coeff_cache.resize(system.atom_num);
	force_coeff_cache.resize(system.atom_num);
	angular_momentum_cache.resize(system.atom_num);
	rot_viscosity_cache.resize(system.atom_num);
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

void SwimForce3d::init_mpi(int thread_count) {
	printf("Swim force: ");
	using_thread = thread_count > 0;
	printf(using_thread ? "Using thread pool\n" : "Using main thread\n");
}

void SwimForce3d::update_ahead(State& state, std::vector<Vec>& F) {

}

void SwimForce3d::update(const State& state, std::vector<Vec>& force_buffer) {
	for (int i = 0; i < force_buffer.size(); ++i) {
		if (!group_cache[i]) continue;
		angular_momentum_cache[i][0] = (rand_uniform() - 0.5) * torque_coeff_cache[i];
		angular_momentum_cache[i][1] = (rand_uniform() - 0.5) * torque_coeff_cache[i];
		angular_momentum_cache[i][2] = (rand_uniform() - 0.5) * torque_coeff_cache[i];
	}

	for (int i = 0; i < force_buffer.size(); i++) {
		direct_cache[i] += direct_cache[i].cross(angular_momentum_cache[i] * rot_coeff_cache[i]);
        direct_cache[i] /= sqrt(direct_cache[i].norm2());
	}

	for (int i = 0; i < force_buffer.size(); i++) {
		if (!group_cache[i])continue;
		force_buffer[i] = direct_cache[i] * force_coeff_cache[i];
	}

}

void SwimForce3d::mp_update(FixedThreadPool& pool, const State& state, std::vector<Vec>& force_buffer) {
	if (using_thread) {
		pool.submit(std::bind(&SwimForce3d::update, this, std::ref(state), std::ref(force_buffer)));
	}
	else {
		update(state, force_buffer);
	}
}

void SwimForce3d::update_cache(const System& system, const Context& context) {

	// group cache
	for (int i = 0; i < system.atom_num; i++) {
		if (system.atom_type[i] == 1)group_cache[i] = 1;
		else group_cache[i] = 0;
	}

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

double SwimForce3d::compute_pressure(const State& state, const std::vector<Vec>& force_buffer) {
	return 0.0;
}

double SwimForce3d::compute_energy(const State& state) {
	return 0.0;
}

#endif // THREE_DIMENSION