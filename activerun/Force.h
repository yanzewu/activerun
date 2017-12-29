#pragma once

#ifndef FORCE_H
#define FORCE_H



#include <math.h>
#include <algorithm>
#include <numeric>

#include "threadpool.h"
#include "dict.h"
#include "mathutil.h"
#include "resources.h"

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

    int* group_cache;

	BrownianForce();

	~BrownianForce();

	void init(const Dict&, System&);

	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec>& force_buffer);

	void update(const State&, std::vector<Vec>& force_buffer);

	void mp_update(FixedThreadPool&, const State& state, std::vector<Vec>& force_buffer);

	void update_cache(const System&, const Context&);

    double compute_temperature(const std::vector<Vec>& force)const;

private:

    std::vector<double> force_coeff_cache;
    std::vector<Vec> random_cache;

};



class SwimForce : public Force {
public:

    double temperature;
    int my_type;

	bool using_thread;
    bool brownian_rotation;

    int* group_cache;


	SwimForce();

	void init(const Dict&, System& system);

	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec>& force_buffer);

	void update(const State&, std::vector<Vec>& force_buffer);

	void mp_update(FixedThreadPool&, const State&, std::vector<Vec>& foce_buffer);

	void update_cache(const System&, const Context&);

private:

    std::vector<double> Pe_R;
	std::vector<double> Pe_S;

	std::vector<double> angle_cache;            // angle of atom
    std::vector<double> angular_momentum_cache; // angular momentum of atom
    std::vector<double> force_coeff_cache;      // swim force
    std::vector<double> rot_coeff_cache;        // dt/zetaR
    std::vector<double> torque_coeff_cache;     // magnitude of torque. T=T0*rand(-0.5,0.5)

};


class SwimForce3d : public Force {
public:

	double temperature;
	int my_type;

	bool using_thread;
	bool brownian_rotation;

	int* group_cache;

	SwimForce3d();

	void init(const Dict&, System& system);

	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec3>& force_buffer);

	void update(const State&, std::vector<Vec3>& force_buffer);

	void mp_update(FixedThreadPool&, const State&, std::vector<Vec3>& foce_buffer);

	void update_cache(const System&, const Context&);

private:

	std::vector<double> Pe_R;
	std::vector<double> Pe_S;

	std::vector<Vec3> angular_momentum_cache;   // angular momentum
	std::vector<Vec3> direct_cache;             // direction of axis
	std::vector<double> force_coeff_cache;      // swim force
    std::vector<double> rot_coeff_cache;        // dt/zetaR
	std::vector<double> torque_coeff_cache;     // magnitude of torque

};


class MorseForce : public Force {
public:

#ifdef PRESSURE_BREAKDOWN
    std::vector<double> pressure;
    std::vector<double> energy;

    std::vector<std::vector<double> > pressure_cache;
    std::vector<std::vector<double> > energy_cache;

#endif

	/* Fixed Data */


	int pool_size;
    double cutoff_relative;
    double cutoff_global;
    double kappa;
    double Um;

	/* Cached Data */

    int* group_cache;

	MorseForce();

	void init(const Dict&, const System&);

	// initialize parallel components. must be called after init.
	void init_mpi(int thread_count);

    // initialize cache for mp
	void update_ahead(double compute_pe);

	// update, submit job in parallel
	void mp_update(FixedThreadPool&, const State&, const Context&);

	// update in main thread
	void update(const State&, const Context&);

	// collect and wait job in parallel
	void update_later(std::vector<Vec>& F);

	void update_cache(const System& system, const Context&);
    
    double max_cutoff()const;

    double potential_energy()const;

private:
	// calculate a btch of column: INSIDE MULTITHREAD
	void update_batch(int start, int end, const State*, const PBCInfo*, const NeighbourList*, std::vector<size_t>* neigh_cache, std::vector<Vec>* force_cache, int tid);

	// calculate single column: INSIDE MULTITHREAD
	void update_column(int i, const State*, const PBCInfo*, const NeighbourList*, std::vector<size_t>* neigh_cache, std::vector<Vec>* force_cache, int tid);

	// calculate a pair
	void update_pair(size_t id1, size_t id2, const Vec& d, const State&, std::vector<Vec>& force_cache, int tid);

	double pair_force_div_r(double r, double exp_cache);

	double pair_energy(double r, double exp_cache);

	double cutoff_global2;	// square of cutoff_global

	std::vector<double> atom_radius_cache;

	std::vector<std::vector<size_t> > neigh_cache;
	std::vector<std::vector<Vec> > force_cache;
	std::vector<std::vector<Vec2d> > pair_cache;

    std::vector<double> pe_cache;

	bool calculate_energy;
};


#endif // !FORCE_H
