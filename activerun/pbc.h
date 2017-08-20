#pragma once

#ifndef SIMUTIL_PBC_H
#define SIMUTIL_PBC_H

#include <stdio.h>
#include <vector>
#include "vec.h"

#include "dimension.h"


class PBCInfo {
public:

	std::vector<Vec> image_cache;

	bool compute_image;

	void init(size_t atom_num, const Vec& box, bool compute_image) {
		box_cache = box;
		half_box_cache = box_cache * 0.5;
		if (this->compute_image = compute_image) {
			location_cache.resize(atom_num);
			image_cache.resize(atom_num);
		}
	}

	// write result into image
	void update(std::vector<Vec>& pos) {
		if (!compute_image) {
			for (auto& p : pos) {
				for (int i = 0; i < DIMENSION; i++) {
					if (p[i] >= box_cache[i]) p[i] -= box_cache[i];
					else if (p[i] < 0) p[i] += box_cache[i];

				}
			}
		}
		else {
			auto p = pos.begin();
			auto loc = location_cache.begin();

			for (; p != pos.end(); p++, loc++) {
				for (int i = 0; i < DIMENSION; i++){
					if ((*p)[i] >= box_cache[i]) { (*p)[i] -= box_cache[i]; (*loc)[i]++; }
					else if ((*p)[i] < 0) { (*p)[i] += box_cache[i]; (*loc)[i] --; }
				}
			}
		}
		

		// most important checkpoint
		for (auto& p : pos) {
			if (p[0] >= box_cache[0] || p[0] < 0
				|| p[1] >= box_cache[1] || p[1] < 0 
#ifdef THREE_DIMENSION
				|| p[2] >= box_cache[2] || p[2] < 0
#endif // THREE_DIMENSION
				) {
				
				fprintf(stderr, "Atom %d missing at %g, %g\n", (int)(&p - &pos[0]), p[0], p[1]);
				throw std::out_of_range("Atom missing");

			}
		}

	}

	// correct the pair distance by pbc
	inline void wrap_pair(Vec& v)const {
		v[0] = v[0] > half_box_cache[0] ? v[0] - box_cache[0] : v[0] < -half_box_cache[0] ? v[0] + box_cache[0] : v[0];
		v[1] = v[1] > half_box_cache[1] ? v[1] - box_cache[1] : v[1] < -half_box_cache[1] ? v[1] + box_cache[1] : v[1];
#ifdef THREE_DIMENSION
		v[2] = v[2] > half_box_cache[2] ? v[2] - box_cache[2] : v[2] < -half_box_cache[2] ? v[2] + box_cache[2] : v[2];

#endif // 3D

	}

	// update image cache
	void update_image(const std::vector<Vec>& pos) {
		for (size_t i = 0; i < pos.size(); i++) {
			image_cache[i] = pos[i] + location_cache[i].to_float() * box_cache;
		}
	}

    // reset location to main box
    void reset_location() {
        location_cache.assign(location_cache.size(), Vecd());
    }

private:
	Vec box_cache;
	Vec half_box_cache;
	std::vector<Vecd> location_cache;


};


#endif // SIMUTIL_PBC_H

