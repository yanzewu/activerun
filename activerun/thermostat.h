#pragma once
#include "System.h"

class ThermoStat {
public:

	std::vector<Vec> init_pos;
	std::vector<double> pressure_cache;
	double temperature_cache;
	
	std::vector<char> compute_pressure;
	bool compute_temperature;
	double volume;

	void init(const Context& context, const State& state_init) {
		pressure_cache.resize(context.force_buffer.size());
		compute_pressure.resize(context.force_buffer.size());
		
		init_pos = state_init.pos;
	}

	void update(const Context& context, const std::vector<Vec> velocity) {
		for (size_t i = 0; i < context.force_buffer.size(); i++) {
			if (compute_pressure[i]) {
				pressure_cache[i] = 0;
				for (size_t j = 0; j < context.force_buffer[i].size(); j++) {
					pressure_cache[i] += context.force_buffer[i][j].dot(context.pbc.image_cache[j] - init_pos[j]);
				}
				pressure_cache[i] /= (volume * DIMENSION); // 2 here is dimension!
			}
		}

		if (compute_temperature) {
			temperature_cache = 0.0;
			for (const auto& v : velocity) {
				temperature_cache += v.norm2();
			}
			temperature_cache /= (velocity.size() * DIMENSION); // 2 here is dimension!
		}
	}

	void update_cache(const System& system) {
		volume = system.box[0] * system.box[1];
	}
};