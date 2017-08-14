#include "Integrator.h"

LangevinIntegrator::LangevinIntegrator() :
	compute_temperature(false),
	cur_force_buffer(force_buffer[0]),
	last_force_buffer(force_buffer[1])
{

}



void LangevinIntegrator::init(const System& system, const Context& context)

{
	inv_viscosity_cache = system.get_attr("zeta");
	for (auto& fb : force_buffer) {
		fb.resize(system.atom_num);
	}

	velocity_cache.resize(system.atom_num);
}

void LangevinIntegrator::update(State& state, const Context& context) {
	memset(&cur_force_buffer[0], 0, sizeof(Vec)*cur_force_buffer.size());
	for (const auto& fb : context.force_buffer) {
		for (size_t i = 0; i < fb.size(); i++) {
			cur_force_buffer[i] += fb[i];
		}
	}

	if (!compute_temperature) {
		for (int i = 0; i < state.pos.size(); i++) {
			state.pos[i] += (cur_force_buffer[i] * inv_viscosity_cache[i]);
		}
	}
	else {
		for (int i = 0; i < state.pos.size(); i++) {
			velocity_cache[i] = cur_force_buffer[i] * inv_viscosity_cache[i];
		}
		for (int i = 0; i < state.pos.size(); i++) {
			state.pos[i] += velocity_cache[i] * timestep_cache;
		}
	}
}

void LangevinIntegrator::update_first_half(State& state) {
	if (!compute_temperature) {
		for (int i = 0; i < state.pos.size(); i++) {
			state.pos[i] += (last_force_buffer[i] * inv_viscosity_cache[i]) * 0.5;
		}
	}
	else {
		for (int i = 0; i < state.pos.size(); i++) {
			velocity_cache[i] = last_force_buffer[i] * inv_viscosity_cache[i];
		}
		for (int i = 0; i < state.pos.size(); i++) {
			state.pos[i] += velocity_cache[i] * timestep_cache * 0.5;
		}
	}
}

void LangevinIntegrator::update_last_half(State& state, const Context& context) {
	memset(&cur_force_buffer[0], 0, sizeof(Vec)*cur_force_buffer.size());
	for (const auto& fb : context.force_buffer) {
		for (size_t i = 0; i < fb.size(); i++) {
			cur_force_buffer[i] += fb[i];
		}
	}

	if (!compute_temperature) {
		for (int i = 0; i < state.pos.size(); i++) {
			state.pos[i] += (cur_force_buffer[i] * inv_viscosity_cache[i]) * 0.5;
		}
	}
	else {
		for (int i = 0; i < state.pos.size(); i++) {
			velocity_cache[i] = cur_force_buffer[i] * inv_viscosity_cache[i];
		}
		for (int i = 0; i < state.pos.size(); i++) {
			state.pos[i] += velocity_cache[i] * timestep_cache * 0.5;
		}
	}

	auto& tmp = cur_force_buffer;
	cur_force_buffer = last_force_buffer;
	last_force_buffer = tmp;

}

double LangevinIntegrator::update_temperature() {
	double temperature = 0.0;
	for (size_t i = 0; i < velocity_cache.size(); i++) {
		temperature += velocity_cache[i].norm2() / inv_viscosity_cache[i];
	}
	return temperature * 0.25 * timestep_cache / velocity_cache.size(); // 0.25-2d!
}

void LangevinIntegrator::update_cache(const Context& context) {
	timestep_cache = context.timestep;
	if (!compute_temperature) {
		for (auto& ivc : inv_viscosity_cache) {
			ivc = context.timestep / ivc;
		}
	}
	else {
		for (auto& ivc : inv_viscosity_cache) {
			ivc = 1.0 / ivc;
		}
		timestep_cache = context.timestep;
	}

}
