#pragma once

#ifndef ACTIVERUN_THERMOSTAT_H
#define ACTIVERUN_THERMOSTAT_H

#include "System.h"
#include "Context.h"
#include "Serializer.h"


    // Manage thermostat 
class ThermoCounter {
public:

    ThermoCounter() {

    }

    void init(const json& thermo_param) {
        thermo_step = thermo_param["thermo_step"];          // step of thermo output
        thermo_start = thermo_param["thermo_start"];        // start of thermo to avoid inequlibrium
        thermo_sample_range = thermo_param["average_range"];// sample range around output time
        thermo_sample_step = thermo_param["average_step"];  // sample step
        thermo_sample_count = thermo_sample_range % thermo_sample_step == 0 ?
            thermo_sample_range / thermo_sample_step + 1 :
            thermo_sample_range / thermo_sample_step;       // sample number
    }

    void check_start(size_t step_begin){
        if (thermo_start < step_begin) {
            printf("Warning: thermo_start too small (%zd), using first step (%zd) instead.\n", thermo_start, step_begin);
            thermo_start = step_begin;
        }
    }

    int sample_count()const {
        return thermo_sample_count;
    }

    bool is_thermo_init_step(size_t step)const {
        return step == thermo_start;
    }

    bool is_thermo_sample_step(size_t step)const {
        if (step >= thermo_start + (thermo_sample_range + 1) / 2) {
            int step_overhead = (step + (thermo_sample_range + 1) / 2 - thermo_start) % thermo_step;
            return (step_overhead <= thermo_sample_range) & (step_overhead % thermo_sample_step == 0);
        }
        else {
            return false;
        }
    }

    bool is_thermo_output_step(size_t step)const {
        if (step >= thermo_start + (thermo_sample_range + 1) / 2) {
            int step_overhead = (step + (thermo_sample_range + 1) / 2 - thermo_start) % thermo_step;
            return (step_overhead == thermo_sample_range);
        }
        else {
            return false;
        }
    }

    size_t last_thermo_step(size_t step)const {
        return (step / thermo_step) * thermo_step;
    }

private:
    size_t thermo_step;
    size_t thermo_start;
    int thermo_sample_range;
    int thermo_sample_step;
    int thermo_sample_count;
};

#endif