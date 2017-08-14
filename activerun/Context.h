#pragma once

#ifndef ACTIVERUN_CONTEXT_H
#define ACTIVERUN_CONTEXT_H

#include "System.h"
#include "neighlist.h"
#include "pbc.h"

struct Context {
	double timestep;
	size_t total_steps;
	size_t current_step;

	NeighbourList* neigh_list;
	PBCInfo pbc;

	std::vector<std::vector<Vec> > force_buffer;

	Context() : neigh_list(nullptr) {

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

	~Context() {
		delete neigh_list;
	}
};


#endif // !ACTIVERUN_CONTEXT_H
