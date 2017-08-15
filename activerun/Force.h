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

    bool update_position_by_self;
    bool calculate_pressure_by_self;
    bool calculate_pressure;
    bool calculate_energy;

    Force(bool is_direct, bool is_paired) : is_direct(is_direct), is_paired(is_paired) {

    }
};

// Brownian Force
class BrownianForce : public Force {
public:

    int my_type;
    double temperature;

	bool using_thread;

    std::vector<int> group_cache;
    std::vector<double> force_coeff_cache;
    std::vector<Vec> random_cache;
    // compute_pressure_by_self = false
    // integrate_by_self = false
    // compute_energy_by_self = true
    // use_celllist = false

	BrownianForce();

	~BrownianForce();

	void init(const Dict&, System&);

	void update_ahead(State&, std::vector<Vec>& force_buffer);

	void update(const State&, std::vector<Vec>& force_buffer);

	double compute_pressure(const State&, const std::vector<Vec>& force);

	double compute_energy(const State&);

	void update_cache(const System&, const Context&);
};



class SwimForce : public Force {
public:

    double temperature;
    int my_type;

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

    SwimForce() : Force(false, false)
    {

    }

    void init(const InputParameter& input_params, System& system) {

        temperature = input_params.kT;
		group_cache.resize(system.atom_num);

		angle_cache.resize(system.atom_num);
		for (auto& angle : angle_cache) {
			angle = fRand(0, M_PI);
		}

		torque_coeff_cache.resize(system.atom_num);
		force_coeff_cache.resize(system.atom_num);
		angular_momentum_cache.resize(system.atom_num);
		rot_viscosity_cache.resize(system.atom_num);
		rot_coeff_cache.resize(system.atom_num);

        // zeta = 6 pi eta r^2
        if (system.attribute_names.find("zeta") == system.attribute_names.end()) {
            auto& zeta = system.add_attr("zeta");
            system.set_name("size", 2);

            for (size_t i = 0; i < system.atom_num; i++) {
                zeta[i] = 6 * M_PI * input_params.viscosity * system.atom_attributes[2][i];
            }
        }

        brownian_rotation = input_params.brownrot;
        Pe_R.resize(system.atom_num, input_params.PeR);
        Pe_S.resize(system.atom_num, input_params.PeS);
    }

    void update_ahead(State& state, Vec2* F) {

    }

    void update(const State& state, std::vector<Vec2>& force_buffer) {
        for (int i = 0; i < force_buffer.size(); ++i) {
            if (!group_cache[i]) continue;
            angular_momentum_cache[i] = torque_coeff_cache[i] * (fRand(-0.5, 0.5));
        }

        for (int i = 0; i < force_buffer.size(); i++) {
            angle_cache[i] += angular_momentum_cache[i] * rot_coeff_cache[i];
        }

        for (int i = 0; i < force_buffer.size(); i++) {
            force_buffer[i] = Vec2(cos(angle_cache[i]), sin(angle_cache[i])) * force_coeff_cache[i];
        }

    }

    void update_cache(const System& system, const Context& context) {

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
                24.0 * temperature * zeta_R[i] * zeta_R[i] / tau_R[i] / context.timestep);
        }
    }

    double compute_pressure(const State& state, const Context& context, const std::vector<Vec2>& force_buffer) {
        double pressure;
        for (int i = 0; i < force_buffer.size(); i++) {
            if (!group_cache[i])continue;
            pressure += force_buffer[i].dot(context.pbc.image_cache[i] - image_begin[i]);
        }
        return pressure;
    }

    double compute_energy(const State& state) {

    }
};



class MorseForce : public Force {
public:

	/* Fixed Data */

    int my_type;

	int pool_size;
    double cutoff_relative;
    double cutoff_global;
    double kappa;
    double Um;

	/* Cached Data */

    double energy_cache;
	double cutoff_global2;	// square of cutoff_global

    std::vector<int> group_cache;
    std::vector<double> atom_radius_cache;

	std::vector<std::vector<size_t> > neigh_cache;
	std::vector<std::vector<Vec> > force_cache;

	MorseForce();

	void init(const Dict&, const System&);

	// initialize parallel components. must called after init.
	void init_mpi(int thread_count);

	void update_ahead(State&, std::vector<Vec>&);

	// update, submit job in parallel
	void update(FixedThreadPool&, const State&, const PBCInfo&, const NeighbourList&);

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


};


#endif // !FORCE_H
