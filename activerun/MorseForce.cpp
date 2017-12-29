
#include "Force.h"
#include "arrayutil.h"
#include "vecutil.h"


MorseForce::MorseForce() {
	this->is_direct = false;
	this->is_paired = true;
	this->is_potential_force = true;
}

void MorseForce::init(const Dict& params, const System& system) {

	logger->write_all("\nInitializing Morse potential\n\n");

	try {
		kappa = params.at("kappa");
		Um = params.at("Um");
	}
	catch (const std::out_of_range&) {
		fprintf(stderr, "Error: kappa or Um not found\n");
		throw;
	}

	cutoff_relative = params.get("cutoff", 1.25);
    cutoff_global = vec_max(system.get_attr("size")) * cutoff_relative;

	logger->write_all("kappa=%.5g\nUm=%.5g\ncutoff=%.4f\n", kappa, Um, cutoff_relative);

#ifdef PRESSURE_BREAKDOWN
    pressure.resize(3);
    energy.resize(3);
#endif 
}

void MorseForce::init_mpi(int thread_count) {

	logger->write_all("Morse potential: ");
	pool_size = thread_count;
    int cache_size;
	if (pool_size > 0) {
		logger->write_all("Using %d threads\n", pool_size);
        cache_size = pool_size;
	}
	else {
		logger->write_all("Using main thread\n");
        cache_size = 1;
	}

	force_cache.resize(cache_size);
	neigh_cache.resize(cache_size);
    pe_cache.resize(cache_size);

	for (auto& nc : neigh_cache) {
		nc.reserve(64);
	}

#ifdef PRESSURE_BREAKDOWN
    pressure_cache.resize(cache_size);
    energy_cache.resize(cache_size);

    for (auto& pc : pressure_cache)pc.resize(3);
    for (auto& ec : energy_cache)ec.resize(3);
#endif

}

void MorseForce::update_ahead(double compute_pe) {

	for (auto& fc : force_cache) {
        std::fill(fc.begin(), fc.end(), Vec());
	}
    this->calculate_energy = compute_pe;
    if (compute_pe) {
#ifdef PRESSURE_BREAKDOWN
        for (auto& ec : energy_cache) vec_reset(ec);
        for (auto& pc : pressure_cache) vec_reset(pc);
        vec_reset(pressure);
        vec_reset(energy);
#else
        vec_reset(pe_cache);
#endif 
    }
}

void MorseForce::mp_update(FixedThreadPool& pool, const State& state, const Context& context) {


	if (pool_size >= 1) { // using parallel

		for (int i = 0; i < pool_size; i++) {
			pool.submit(std::bind(&MorseForce::update_batch,
				this,
				1 + context.neigh_list->box_num[0] * i / pool_size, 1 + context.neigh_list->box_num[0] * (i + 1) / pool_size,
				&state, &context.pbc, context.neigh_list, &neigh_cache[i], &force_cache[i], i), true);
		}
	}
	else {
		update(state, context);
	}
}

void MorseForce::update(const State& state, const Context& context) {
	update_batch(1, context.neigh_list->box_num[0] + 1,
		&state, &context.pbc, context.neigh_list, &neigh_cache[0], &force_cache[0], 0);
}

void MorseForce::update_later(std::vector<Vec>& F) {
	for (const auto& fc : force_cache) {
		array_add(fc, F);
	}

#ifdef PRESSURE_BREAKDOWN
    if (calculate_energy) {
        for (const auto& pc : pressure_cache) array_add_double(&pc[0], &pressure[0], 3);
        for (const auto& ec : energy_cache) array_add_double(&ec[0], &energy[0], 3);
    }
#endif 
}

void MorseForce::update_batch(int start, int end, const State* state, const PBCInfo* pbc, const NeighbourList* neigh_list,
	std::vector<size_t>* neigh_cache, std::vector<Vec>* F, int tid) {
	for (int i = start; i < end; i++) {
		update_column(i, state, pbc, neigh_list, neigh_cache, F, tid);
	}
}

void MorseForce::update_column(int i, const State* state, const PBCInfo* pbc, const NeighbourList* neigh_list, std::vector<size_t>* neigh_cache, std::vector<Vec>* F, int tid) {

	for (int j = 1; j <= neigh_list->box_num[1]; j++) 
#ifdef THREE_DIMENSION
		for(int k = 1; k <= neigh_list->box_num[2]; k++){
		auto cell1 = &neigh_list->at(i, j, k);
#else
		{
		auto cell1 = &neigh_list->at(i, j);
#endif // THREE_DIMENSION
		neigh_cache->clear();

#ifdef THREE_DIMENSION
		{auto cell2 = &neigh_list->at(i + 1, j - 1, k); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		{auto cell2 = &neigh_list->at(i + 1, j, k); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		{auto cell2 = &neigh_list->at(i + 1, j + 1, k); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		{auto cell2 = &neigh_list->at(i, j + 1, k); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		for(int p = -1; p <= 1; p++)
		for (int q = -1; q <= 1; q++)
			{
				auto cell2 = &neigh_list->at(i + p, j + q, k + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); 

			}

#else
		{auto cell2 = &neigh_list->at(i + 1, j - 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		{auto cell2 = &neigh_list->at(i + 1, j); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		{auto cell2 = &neigh_list->at(i + 1, j + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }
		{auto cell2 = &neigh_list->at(i, j + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end()); }

#endif // THREE_DIMENSION

		
		for (const auto& id1 : *cell1)
			for (const auto& id2 : *neigh_cache) {
				//		if (id2 >= id1) continue;
				auto d = state->pos[id2] - state->pos[id1];
				pbc->wrap_pair(d);
				if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global
#ifdef THREE_DIMENSION
					|| abs(d[2]) > cutoff_global
#endif // THREE_DIMENSION

					)continue;
				update_pair(id1, id2, d, *state, *F, tid);
			}

		for (const auto& id1 : *cell1)
			for (const auto& id2 : *cell1) {
				if (id2 >= id1) continue;
				auto d = state->pos[id2] - state->pos[id1];
				pbc->wrap_pair(d);
				if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global
#ifdef THREE_DIMENSION
					|| abs(d[2]) > cutoff_global
#endif // THREE_DIMENSION
					)continue;
				update_pair(id1, id2, d, *state, *F, tid);
			}
	}

}

void MorseForce::update_pair(size_t id1, size_t id2, const Vec& d, const State& state, std::vector<Vec>& F, int tid) {
	
	double r2 = d.norm2();
	if (r2 > cutoff_global2)return;

	double sum_radius = atom_radius_cache[id1] + atom_radius_cache[id2];

	double cutoff = sum_radius * cutoff_relative;
	if (r2 > cutoff * cutoff)return; //cutoff

	double r = sqrt(r2);


	double exp_cache = exp(-kappa * (r - sum_radius));
#ifdef DEBUG_FORCE
	if (abs(exp_cache) > 1000) {
	    fprintf(stderr, "Warning: Force too large for %d-%d, %.5g\n", id1, id2, exp_cache);
	    fprintf(stderr, "Distance is %.4f, %.4f (%.4f)\n", d[0], d[1], r);
	    fprintf(stderr, "Position is %.4f,%.4f and %.4f,%.4f\n", state.pos[id1][0], state.pos[id1][1], state.pos[id2][0], state.pos[id2][1]);
	    fprintf(stderr, "Radius is %.4f and %.4f\n", atom_radius_cache[id1], atom_radius_cache[id2]);
	    throw std::out_of_range("Force too large");
    }
#endif
	Vec f_temp = d * pair_force_div_r(r, exp_cache) / r;
	F[id1] -= f_temp;
	F[id2] += f_temp;


    if (calculate_energy) {
#ifdef PRESSURE_BREAKDOWN
        int output_id;

        if (group_cache[id1] && group_cache[id2])output_id = 1;
        else if (group_cache[id1] || group_cache[id2])output_id = 2;
        else output_id = 0;

        energy_cache[tid][output_id] += pair_energy(r, exp_cache);
        pressure_cache[tid][output_id] += f_temp.dot(d);
#else
        pe_cache[tid] += pair_energy(r, exp_cache);
#endif
    }
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

double MorseForce::max_cutoff()const {
    return cutoff_global;
}

double MorseForce::potential_energy()const {
    return vec_sum(pe_cache);
}