#pragma once

#include "System.h"
#include "Context.h"

class FixPressureComputer {
public:

    int* group_cache;

    void init() {
        group_cache = nullptr;
    }

    double compute_pressure(const Context& context, int force_id) {
        double pressure = 0.0;

        for (size_t j = 0; j < context.force_buffer[force_id].size(); j++) {
            if (group_cache && !group_cache[j])continue;
            pressure += context.force_buffer[force_id][j].dot(context.pbc.image_cache[j] - init_pos[j]);
        }

        return pressure / (volume * DIMENSION);
    }

    void update_cache(const System& system, const Context& context) {
#ifndef THREE_DIMENSION
        volume = system.box[0] * system.box[1];
#else
        volume = system.box[0] * system.box[1] * system.box[2];
#endif // !THREE_DIMENSION

        init_pos = context.pbc.image_cache;
    }

private:
    std::vector<Vec> init_pos;
    double volume;
};

class KineticEnergyComputer {
public:

    void init() {

    }

    double compute_temperature(const std::vector<Vec>& velocity) {
        double temperature = 0.0;
        for (const auto& v : velocity) {
            temperature += v.norm2();
        }
        return temperature / (velocity.size() * DIMENSION);

    }
};
