#pragma once

#ifndef FORCE_H
#define FORCE_H



#include <math.h>
#include <algorithm>
#include <numeric>

#include "threadpool.h"
#include "dict.h"
#include "mathutil.h"

#include "System.h"
#include "Context.h"

class Force {
public:
    bool is_direct; // applied directly, without integrator
    bool is_paired; // applied pair-wise, may need celllist
    bool is_potential_force;	// calculate potential energy

    Force() {

    }
};

// Brownian Force
class BrownianForce : public Force {
public:

    double temperature;

	bool using_thread;

    std::vector<int> group_cache;
    std::vector<double> force_coeff_cache;
    std::vector<Vec> random_cache;

	BrownianForce();

	~BrownianForce();

	void init(const Dict&, System&);

	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec>& force_buffer);

	void update(const State&, std::vector<Vec>& force_buffer);

	void mp_update(FixedThreadPool&, const State& state, std::vector<Vec>& force_buffer);

	double compute_pressure(const State&, const std::vector<Vec>& force);

	double compute_energy(const State&);

	void update_cache(const System&, const Context&);
};



class SwimForce : public Force {
public:

    double temperature;
    int my_type;

	bool using_thread;

    bool brownian_rotation;

    std::vector<double> Pe_R;
	std::vector<double> Pe_S;

    Vec2* image_begin;

    std::vector<int> group_cache;

	std::vector<double> angle_cache;

    std::vector<double> torque_coeff_cache;
    std::vector<double> angular_momentum_cache;
    std::vector<double> force_coeff_cache;
    std::vector<double> rot_viscosity_cache;
    std::vector<double> rot_coeff_cache;

	SwimForce();

	void init(const Dict&, System& system);

	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec>& force_buffer);

	void update(const State&, std::vector<Vec>& force_buffer);

	void mp_update(FixedThreadPool&, const State&, std::vector<Vec>& foce_buffer);

	void update_cache(const System&, const Context&);

	double compute_pressure(const State&, const std::vector<Vec>& force_buffer);

	double compute_energy(const State&);
};



class MorseForce : public Force {
public:

	/* Fixed Data */

	bool calculate_energy;

	int pool_size;
    double cutoff_relative;
    double cutoff_global;
    double kappa;
    double Um;

	/* Cached Data */

    std::vector<int> group_cache;

	MorseForce();

	void init(const Dict&, const System&);

	// initialize parallel components. must called after init.
	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec>&);

	// update, submit job in parallel
	void mp_update(FixedThreadPool&, const State&, const PBCInfo&, const NeighbourList&);

	// update in main thread
	void update(const State&, const PBCInfo&, const NeighbourList&);

	// collect and wait job in parallel
	void update_later(std::vector<Vec>& F);

	void update_cache(const System& system, const Context& context);

	double compute_pressure(const State& state);

	double compute_energy(const State& state);

private:
	// calculate a btch of column: INSIDE MULTITHREAD
	void update_batch(int start, int end, const State*, const PBCInfo*, const NeighbourList*, std::vector<size_t>* neigh_cache, std::vector<Vec>* force_cache);

	// calculate single column: INSIDE MULTITHREAD
	void update_column(int i, const State*, const PBCInfo*, const NeighbourList*, std::vector<size_t>* neigh_cache, std::vector<Vec>* force_cache);

	// calculate a pair
	void update_pair(size_t id1, size_t id2, const Vec& d, const State&, std::vector<Vec>& force_cache);

	double pair_force_div_r(double r, double exp_cache);

	double pair_energy(double r, double exp_cache);


	double energy_cache;
	double cutoff_global2;	// square of cutoff_global

	std::vector<double> atom_radius_cache;

	std::vector<std::vector<size_t> > neigh_cache;
	std::vector<std::vector<Vec> > force_cache;
	std::vector<std::vector<Vec2d> > pair_cache;
};


#endif // !FORCE_H
