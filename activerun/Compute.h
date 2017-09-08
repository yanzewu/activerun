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
        volume = system.volume();
        init_pos = context.pbc.image_cache;
    }

private:
    std::vector<Vec> init_pos;
    double volume;
};

class PressureComputer {
public:

    bool using_pbc;
    int* group_cache;

    void init(double using_pbc) {
        this->using_pbc = using_pbc;
        group_cache = nullptr;
    }

    double compute_pressure(const Context& context, const State& state, int force_id) {

        double pressure = 0.0;

        for (size_t j = 0; j < context.force_buffer[force_id].size(); j++) {
            if (group_cache && !group_cache[j])continue;
            pressure += context.force_buffer[force_id][j].dot(state.pos[j]);
        }

        std::vector<size_t> ghost;
        std::vector<Vecd> offset;

        context.neigh_list->compute_ghost_pos(ghost, offset);
        for (size_t i = 0; i < ghost.size(); i++) {
            if (group_cache && !group_cache[ghost[i]])continue;
            pressure += context.force_buffer[force_id][ghost[i]].dot(state.pos[ghost[i]] + offset[i].to_float() * box);
        }

        return pressure / (DIMENSION * volume);
    }

    void update_cache(const System& system, const Context& context) {
        volume = system.volume();
        box = system.box;
    }

private:

    double volume;
    Vec box;
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
