#pragma once

#ifndef FORCE_H
#define FORCE_H



#include <math.h>
#include <algorithm>
#include <numeric>

#include "threadpool.h"

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

	void init(const InputParameter& input_params, System& system);

	void update_ahead(State& state, std::vector<Vec>& force_buffer);

	void update(const State& state, std::vector<Vec>& force_buffer);

	double compute_pressure(const State& state, const std::vector<Vec>& force);

	double compute_energy(const State& state);

	void update_cache(const System& system, const Context& context);
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

    int my_type;

    double cutoff_relative;
    double cutoff_global;
	double cutoff_global2;
    double kappa;
    double Um;

    int* group_cache;
    std::vector<double> atom_radius_cache;

    double energy_cache;

	bool using_thread;
	int pool_size;

	std::vector<std::vector<size_t> > neigh_cache;
	std::vector<std::vector<Vec2> > force_cache;

    MorseForce() : Force(false, true) {
		pool_size = 3;
    }

    void init(const InputParameter& input_params, const System& system) {
        kappa = input_params.kappa;
        Um = input_params.Um;

        cutoff_relative = 1.0 + input_params.Rg * 2;
        cutoff_global = *std::max_element(system.get_attr("size").begin(), system.get_attr("size").end()) * cutoff_relative;
		cutoff_global2 = cutoff_global * cutoff_global;
    }

    void update_ahead(State& state, std::vector<Vec2>& F) {

    }

    void update(FixedThreadPool& pool, const State& state, const PBCInfo& pbc, const NeighbourList& neigh_list, std::vector<Vec2>& F) {
        energy_cache = 0.0;

		for (auto& fc : force_cache) {
			memset(&fc[0], 0, sizeof(Vec)*fc.size());
		}

/*		for (int i = 0; i < pool_size; i++) {
			pool.submit(std::bind(&MorseForce::update_batch_column,
				this,
				1 + neigh_list.box_num[0] * i / pool_size, 1 + neigh_list.box_num[0] * (i+1) / pool_size,
				&state, &pbc, &neigh_list, &neigh_cache[i], &force_cache[i]));
		}*/

		for (int i = 1; i <= neigh_list.box_num[0]; i++) {
			update_column(i, &state, &pbc, &neigh_list, &neigh_cache[0], &F);
//			pool.submit(std::bind(&MorseForce::update_column, this, i, &state, &pbc, &neigh_list, &neigh_cache[i-1], &force_cache[i-1]));

/*			for (int j = 0; j < neigh_list.box_num[1]; j++) {
				for (int nb_i = i - 1; nb_i <= i + 1; nb_i++)
					for (int nb_j = j - 1; nb_j <= j + 1; nb_j++) {
						for (const auto& id1 : neigh_list(i, j))
							for (const auto& id2 : neigh_list(nb_i, nb_j)) {
								if (id1 >= id2)continue;
								auto d = state.pos[id2] - state.pos[id1];
								pbc.wrap_pair(d);
								if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global)continue;
								update_pair(id1, id2, d, state, F);
							}
					}
			}*/
		}
//		pool.wait();
    }

	void update_later(std::vector<Vec2>& F) {
		for (const auto& fc : force_cache) {
			for (size_t i = 0; i < F.size(); i++) {
				F[i] += fc[i];
			}
		}

	}

	void update_batch_column(int start, int end, const State* state, const PBCInfo* pbc, const NeighbourList* neigh_list,
		std::vector<size_t>* neigh_cache, std::vector<Vec2>* F) {
		for (int i = start; i < end; i++) {
			update_column(i, state, pbc, neigh_list, neigh_cache, F);
		}
	}

	void update_column(int i, const State* state, const PBCInfo* pbc, const NeighbourList* neigh_list, std::vector<size_t>* neigh_cache, std::vector<Vec2>* F) {
		
		for (int j = 1; j <= neigh_list->box_num[1]; j++) {
			auto cell1 = &neigh_list->at(i, j);
			neigh_cache->clear();
	/*		for (int nb_i = i - 1; nb_i <= i + 1; nb_i++)
			int nb_i = i + 1;
			for (int nb_j = j - 1; nb_j <= j + 1; nb_j++) {
				auto cell2 = &neigh_list->at(nb_i, nb_j);
				neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
			}*/

			const std::vector<size_t>* cell2;
			cell2 = &neigh_list->at(i + 1, j - 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
			cell2 = &neigh_list->at(i + 1, j    ); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
			cell2 = &neigh_list->at(i + 1, j + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());
			cell2 = &neigh_list->at(i    , j + 1); neigh_cache->insert(neigh_cache->end(), cell2->begin(), cell2->end());

			for (const auto& id1 : *cell1)
			for (const auto& id2 : *neigh_cache) {
		//		if (id2 >= id1) continue;
				auto d = state->pos[id2] - state->pos[id1];
				pbc->wrap_pair(d);
				if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global)continue;
				update_pair(id1, id2, d, *state, *F);

			}	

			for(const auto& id1 : *cell1)
			for (const auto& id2 : *cell1) {
				if (id2 >= id1) continue;
				auto d = state->pos[id2] - state->pos[id1];
				pbc->wrap_pair(d);
				if (abs(d[0]) > cutoff_global || abs(d[1]) > cutoff_global)continue;
				update_pair(id1, id2, d, *state, *F);
			}
		}

	}

    void update_pair(int id1, int id2, const Vec2& d, const State& state, std::vector<Vec2>& F) {

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
        Vec2 f_temp = d * pair_force_div_r(r, exp_cache) / r;
        F[id1] -= f_temp;
        F[id2] += f_temp;

        if (calculate_energy)energy_cache += pair_energy(r, exp_cache);
    }

    double pair_force_div_r(double r, double exp_cache) {
        return (2 * kappa * Um * exp_cache)*(exp_cache - 1.0);	// negative - attract, positive -repel
    }

    double pair_energy(double r, double exp_cache) {
        return Um * exp_cache * (exp_cache - 2.0);
    }

    void update_cache(const System& system, const Context& context) {
        atom_radius_cache = system.get_attr("size");
        for (size_t i = 0; i < atom_radius_cache.size(); i++) {
            atom_radius_cache[i] *= 0.5;
        }
		force_cache.resize(pool_size);
		neigh_cache.resize(pool_size);
		for (auto& fc : force_cache) {
			fc.resize(system.atom_num);
		}
		for (auto& nc : neigh_cache) {
			nc.reserve(64);
		}
    }

    double compute_pressure(const State& state) {

    }

    double compute_energy(const State& state) {
        return energy_cache;
    }
};

/*

class AOForce : public Force {
public:

    int my_type;

    double cutoff_global;
    double potential_depth;
    double Um;

    int* group_cache;

    double energy_cache;

    AOForce() : Force(false, true) {

    }

    void initialize(const InputParameter& input_params, const System& system) {

        Um = input_params.Um;

        potential_depth = input_params.Rg * 2;
        cutoff_global = *std::max_element(system.atom_type.begin(), system.atom_type.end()) * (1.0 + potential_depth);
    }

    void update_ahead(State& state, Vec2* F) {

    }

    void update(const State& state, Vec2* F) {
        energy_cache = 0.0;
    }

    void update_cell(FixArray<int>* cell1, FixArray<int>* cell2, const State& state, Vec2* F) {
        for (const int& id1 : *cell1)
            for (const int& id2 : *cell2) {
                if (id1 < id2) {
                    update_pair(id1, id2, state, F);
                }
            }
    }

    void update_pair(int id1, int id2, const State& state, Vec2* F) {
        Vec2 d(state.x[id2] - state.x[id1], state.y[id2] - state.y[id1]);
        d[0] = fmod(d[0], state.cache.box[0]);
        d[1] = fmod(d[0], state.cache.box[1]);

        if (d[0] > cutoff_global || d[1] > cutoff_global)return;

        double r = d.norm();
        double sum_radius = state.cache.atom_radius[id1] + state.cache.atom_radius[id2];

        if (r > sum_radius + potential_depth)return; //cutoff

        Vec2 f_temp = d * -pair_force_div_r(r);
        F[id1] += f_temp;
        F[id2] -= f_temp;

        if (calculate_energy)energy_cache += pair_energy(r);
    }

    double pair_force_div_r(double r) {
        return 2.0 * Um / (r * potential_depth * potential_depth) * ();
    }

    double pair_energy(double r) {

    }

    void update_cache(const System& system, const Context& context) {

    }

    double compute_pressure(const State& state) {

    }

    double compute_energy(const State& state) {
        return energy_cache;
    }

};

class PotentialFreeForce {
public:
    void update_ahead(State& state, Vec2* F) {

        std::deque<int[2]> queue;

        for (int cell_x = 0; cell_x < state.celllist.cellcount_x; cell_x++)
            for (int cell_y = 0; cell_y < state.celllist.cellcount_y; cell_y++) {
                for (int i = 0; i <= 1; i++)
                    for (int j = 0; j <= 1; j++) {
                        int neighbour_x = (cell_x + i) % state.celllist.cellcount_x;
                        int neighbour_y = (cell_y + i) % state.celllist.cellcount_y;

                        if (update_cell_ahead(&state.celllist.neighbour[cell_x][cell_y],
                            &state.celllist.neighbour[neighbour_x][neighbour_y], state, F
                        )) {
                            queue.push_back({ cell_x, cell_y });
                            queue.push_back({ neighbour_x, neighbour_y });
                        }
                    }
            }

        while (!queue.empty()) {
            int cell_x = queue.front()[0];
            int cell_y = queue.front()[1];

            queue.pop_front();

            for (int i = 0; i <= 1; i++)
                for (int j = 0; j <= 1; j++) {
                    int neighbour_x = (cell_x + i) % state.celllist.cellcount_x;
                    int neighbour_y = (cell_y + i) % state.celllist.cellcount_y;

                    if (update_cell_ahead(&state.celllist.neighbour[cell_x][cell_y],
                        &state.celllist.neighbour[neighbour_x][neighbour_y], state, F
                    )) {
                        queue.push_back({ cell_x, cell_y });
                        queue.push_back({ neighbour_x, neighbour_y });
                    }
                }
        }
    }

    bool update_cell_ahead(FixArray<int>* cell1, FixArray<int>* cell2, State& state, Vec2* F) {

    }
};

*/

#endif // !FORCE_H
