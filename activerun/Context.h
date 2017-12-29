#pragma once

#ifndef ACTIVERUN_CONTEXT_H
#define ACTIVERUN_CONTEXT_H

#include "System.h"
#include "neighlist.h"
#include "pbc.h"
#include "threadpool.h"
#include "resources.h"

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
		total_steps = (size_t)param.at("step");
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
		logger->write_all("\nInitializing celllist\n\n");
		logger->write_all("cutoff=%f\n", cutoff);

		neigh_list = new NeighbourList(box, cutoff, true);
		if (DIMENSION == 2) {
			logger->write_all("Total %d x %d boxes\n", neigh_list->box_num[0], neigh_list->box_num[1]);
			logger->write_all("Box size=%.4f, %.4f", neigh_list->unit[0], neigh_list->unit[1]);
		}
		else {
			logger->write_all("Total %d x %d x %d boxes\n", neigh_list->box_num[0], neigh_list->box_num[1], neigh_list->box_num[2]);
			logger->write_all("Box size=%.4f, %.4f, %.4f", neigh_list->unit[0], neigh_list->unit[1], neigh_list->unit[2]);
		}
		if (neigh_list->using_ghost)logger->write_all("Using ghost\n");
	}

	void init_multicore(int thread_count, int pair_force_idx) {
		logger->write_all("\nInitialize parallel\n\n");

		thread_num.resize(force_buffer.size());
		for (size_t i = 0; i < thread_num.size(); i++) {
			if (i == pair_force_idx)thread_num[i] = thread_count - 1;
			else thread_num[i] = 0;
		}

		if (thread_count > 1) {
			logger->write_all("Create pool with %d threads\n", thread_count - 1);
			pool.init(thread_count - 1);
		}
		else {
			logger->write_all("No pool created\n");
		}
	}

    void clear_buffer() {
        for (auto& fb : force_buffer) {
            memset(&fb[0], 0, sizeof(Vec) * fb.size());
        }

    }

	~Context() {
		delete neigh_list;
	}
};


#endif // !ACTIVERUN_CONTEXT_H
