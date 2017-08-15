
#include "Force.h"

MorseForce::MorseForce() : Force(false, true) {

}

void MorseForce::init(const Dict& params, const System& system) {

	printf("\nInitializing Morse potential\n\n");

	try {
		kappa = params.at("kappa");
		Um = params.at("Um");
	}
	catch (const std::out_of_range&) {
		printf("Error: kappa or Um not found\n");
		throw;
	}

	cutoff_relative = params.get("cutoff", 1.25);

	printf("kappa=%.5g\nUm=%.5g\ncutoff=%.4f\n", kappa, Um, cutoff_relative);

	calculate_energy = (bool)params.get("energy", 0);
	printf(calculate_energy ? "Using energy cache\n" : "No energy cache\n");
}

void MorseForce::init_mpi(int thread_count) {
	pool_size = thread_count;
	if (pool_size > 0) {
		force_cache.resize(pool_size);
		neigh_cache.resize(pool_size);
	}
	else {
		force_cache.resize(1);
		neigh_cache.resize(1);
	}

	for (auto& nc : neigh_cache) {
		nc.reserve(64);
	}

}

void MorseForce::update_ahead(State& state, std::vector<Vec>& F) {

}

void MorseForce::mp_update(FixedThreadPool& pool, const State& state, const PBCInfo& pbc, const NeighbourList& neigh_list) {
	energy_cache = 0.0;

	for (auto& fc : force_cache) {
		memset(&fc[0], 0, sizeof(Vec)*fc.size());
	}

	if (pool_size >= 1) { // using parallel

		for (int i = 0; i < pool_size; i++) {
			pool.submit(std::bind(&MorseForce::update_batch,
				this,
				1 + neigh_list.box_num[0] * i / pool_size, 1 + neigh_list.box_num[0] * (i + 1) / pool_size,
				&state, &pbc, &neigh_list, &neigh_cache[i], &force_cache[i]), true);
		}
	}
	else {
		update(state, pbc, neigh_list);
	}
}

void MorseForce::update(const State& state, const PBCInfo& pbc, const NeighbourList& neigh_list) {
	update_batch(1, neigh_list.box_num[0] + 1,
		&state, &pbc, &neigh_list, &neigh_cache[0], &force_cache[0]);
}

void MorseForce::update_later(std::vector<Vec2>& F) {
	for (const auto& fc : force_cache) {
		auto fc_iter = fc.begin();
		auto F_iter = F.begin();

		for (; fc_iter != fc.end();) {
			*(F_iter++) += *(fc_iter++);
		}
	}

}

void MorseForce::update_batch(int start, int end, const State* state, const PBCInfo* pbc, const NeighbourList* neigh_list,
	std::vector<size_t>* neigh_cache, std::vector<Vec2>* F) {
	for (int i = start; i < end; i++) {
		update_column(i, state, pbc, neigh_list, neigh_cache, F);
	}
}

void MorseForce::update_column(int i, const State* state, const PBCInfo* pbc, const NeighbourList* neigh_list, std::vector<size_t>* neigh_cache, std::vector<Vec2>* F) {

	for (int j = 1; j <= neigh_list->box_num[1]; j++) {
		auto cell1 = &neigh_list->at(i, j);
		neigh_cache->clear();
		/*		for (int nb_i = i - 1; nb_i <= i + 1; nb_i++)
		int nb_i = i + 1;
		for (int nb_j = j - 1; nb_j <= j + 1; nb_j++) {
		auto cell2 = &neigh_list->at(nb_i, nb_j);
		neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
		}*/

		/// Related with dimension

		const std::vector<size_t>* cell2;
		cell2 = &neigh_list->at(i + 1, j - 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
		cell2 = &neigh_list->at(i + 1, j); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
		cell2 = &neigh_list->at(i + 1, j + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
		cell2 = &neigh_list->at(i, j + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());

		for (const auto& id1 : *cell1)
			for (const auto& id2 : *neigh_cache) {
				//		if (id2 >= id1) continue;
				auto d = state->pos[id2] - state->pos[id1];
				pbc->wrap_pair(d);
				if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global)continue;
				update_pair(id1, id2, d, *state, *F);

			}

		for (const auto& id1 : *cell1)
			for (const auto& id2 : *cell1) {
				if (id2 >= id1) continue;
				auto d = state->pos[id2] - state->pos[id1];
				pbc->wrap_pair(d);
				if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global)continue;
				update_pair(id1, id2, d, *state, *F);
			}
	}

}

void MorseForce::update_pair(size_t id1, size_t id2, const Vec& d, const State& state, std::vector<Vec>& F) {
	
	double r2 = d.norm2();
	if (r2 > cutoff_global2)return;

	double sum_radius = atom_radius_cache[id1] + atom_radius_cache[id2];

	double cutoff = sum_radius * cutoff_relative;
	if (r2 > cutoff * cutoff)return; //cutoff

	double r = sqrt(r2);


	double exp_cache = exp(-kappa * (r - sum_radius));
	/*		if (abs(exp_cache) > 1000) {
	fprintf(stderr, "Warning: Force too large for %d-%d, %.5g\n", id1, id2, exp_cache);
	fprintf(stderr, "Distance is %.4f, %.4f (%.4f)\n", d[0], d[1], r);
	fprintf(stderr, "Position is %.4f,%.4f and %.4f,%.4f\n", state.pos[id1][0], state.pos[id1][1], state.pos[id2][0], state.pos[id2][1]);
	fprintf(stderr, "Radius is %.4f and %.4f\n", atom_radius_cache[id1], atom_radius_cache[id2]);
	//			throw std::out_of_range("Force too large");
	}*/
	Vec f_temp = d * pair_force_div_r(r, exp_cache) / r;
	F[id1] -= f_temp;
	F[id2] += f_temp;

	if (calculate_energy)energy_cache += pair_energy(r, exp_cache);
}

double MorseForce::pair_force_div_r(double r, double exp_cache) {
	return (2 * kappa * Um * exp_cache)*(exp_cache - 1.0);	// negative - attract, positive -repel
}

double MorseForce::pair_energy(double r, double exp_cache) {
	return Um * exp_cache * (exp_cache - 2.0);
}

void MorseForce::update_cache(const System& system, const Context& context) {

	cutoff_global = *std::max_element(system.get_attr("size").begin(), system.get_attr("size").end()) * cutoff_relative;
	cutoff_global2 = cutoff_global * cutoff_global;


	atom_radius_cache = system.get_attr("size");
	for (size_t i = 0; i < atom_radius_cache.size(); i++) {
		atom_radius_cache[i] *= 0.5;
	}
	for (auto& fc : force_cache) {
		fc.resize(this->atom_radius_cache.size());
	}
}

double MorseForce::compute_pressure(const State& state) {
	return 0.0;
}

double MorseForce::compute_energy(const State& state) {
	return energy_cache;
}

