#pragma once

#include "System.h"
#include "Context.h"

class FixPressureComputer {
public:

    int* group_cache;

	void init();

	double compute_pressure(const Context& context, int force_id);

	void update_cache(const System& system);

	void init_computing(const Context& context);

	void read_restart(const char* filename, Context& context);

	void write_restart(const char* filename, const Context& context)const;


private:
    double volume;
    std::vector<Vec> init_pos;
};

class PressureComputer {
public:

    bool using_pbc;
    int* group_cache;

	void init(double using_pbc);

	double compute_pressure(const Context& context, const State& state, int force_id);

	void update_cache(const System& system, const Context& context);

private:

    double volume;
    Vec box;
};

class KineticEnergyComputer {
public:

	void init();

	double compute_temperature(const std::vector<Vec>& velocity);
};
