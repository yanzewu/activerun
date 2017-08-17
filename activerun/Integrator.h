#pragma once
#include "System.h"
#include "Context.h"

/*
	Langevin Integrator.
*/
class LangevinIntegrator {
public:

    std::vector<Vec> velocity_cache;		// only cached when compute_temperature is set

    bool compute_temperature;	// whether temperature is computed (using different update algorithm)

	LangevinIntegrator();

	// must be called after construct
	void init(const Dict& dict, const System& system, const Context& context);

	// perform simple update (require smaller timestep)
	void update(State& state, const Context& context);

	// first part of double-update. NOTICE: double-update is incompatible with Brownian force
	void update_first_half(State& state);

	// second part of double-update
	void update_last_half(State& state, const Context& context);

	// compute temperature. Only valid when compute_temperature is set.
	double update_temperature();

	// must be called before first run
	void update_cache(const System& system, const Context& context);

private:
    double timestep_cache;
    std::vector<double> inv_viscosity_cache;	// acceleration coeff

    std::vector<Vec> force_buffer[2];		// overall force buffer
	std::vector<Vec>& cur_force_buffer;	// pointer to force_buffer
	std::vector<Vec>& last_force_buffer;	// pointer to last force_buffer

};