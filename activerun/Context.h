#pragma once

#ifndef ACTIVERUN_CONTEXT_H
#define ACTIVERUN_CONTEXT_H

#include "System.h"
#include "neighlist.h"
#include "pbc.h"
#include "threadpool.h"

struct Context {
	double timestep;
	size_t total_steps;
	size_t current_step;

	NeighbourList* neigh_list;
	PBCInfo pbc;

	std::vector<std::vector<Vec> > force_buffer;	// buffer for each force

	std::vector<int> thread_num;					// threads for each force
	FixedThreadPool pool;

	Context() : neigh_list(nullptr) {

	}

	void init_timestep(const Dict& param) {
		timestep = param.at("timestep");
		total_steps = (size_t)param.at("steps");
	}

	void init_force_buffer(const System& system, int force_num) {
		force_buffer.resize(force_num);
		for (auto& b : force_buffer) {
			b.resize(system.atom_num);
		}
	}

	void init_pbc(const System& system, bool compute_image) {
		pbc.init(system.atom_num, system.box, compute_image);
	}

	void init_neighlist(const Vec& box, double cutoff) {
		neigh_list = new NeighbourList(box, cutoff, true);
	}

	void init_multicore(int thread_count, int pair_force_idx) {
		printf("\nInitialize parallel\n\n");

		thread_num.resize(force_buffer.size());
		for (size_t i = 0; i < thread_num.size(); i++) {
			if (i == pair_force_idx)thread_num[i] = thread_count - 1;
			else thread_num[i] = 0;
		}

		if (thread_count > 1) {
			printf("Create pool with %d threads\n", thread_count - 1);
			pool.init(thread_count - 1);
		}
		else {
			printf("No pool created\n");
		}
	}

	~Context() {
		delete neigh_list;
	}
};


#endif // !ACTIVERUN_CONTEXT_H
