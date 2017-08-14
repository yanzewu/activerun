#pragma once

#ifndef ACTIVERUN_DUMPER
#define ACTIVERUN_DUMPER



#include <stdio.h>
#include <vector>
#include "System.h"

class Dumper {
public:
    FILE* ofile;

    Dumper(const char* filename) {
        ofile = fopen(filename, "w");
        if (!ofile) {
            printf("Cannot open %s", filename);
            // throw;
        }
    }

    void write(const char* str) {
        fprintf(ofile, str);
    }

    ~Dumper() {
        fclose(ofile);
    }
};

class TrajDumper : public Dumper {
public:

    TrajDumper(const char* filename) : Dumper(filename) {

    }


    void dump(const System& system, const State& state, int step) {
        fprintf(ofile, "ITEM: TIMESTEP\n");
        fprintf(ofile, "%d\n", step);
        fprintf(ofile, "ITEM: NUMBER OF ATOMS\n");
        fprintf(ofile, "%i\n", system.atom_num);
        fprintf(ofile, "ITEM: BOX BOUNDS\n");
        fprintf(ofile, "0 %f\n", system.box[0]);
        fprintf(ofile, "0 %f\n", system.box[1]);
        fprintf(ofile, "-0.25 0.25\n");
        fprintf(ofile, "ITEM: ATOMS id mol type x y z\n");
        for (int i = 0; i < system.atom_num; ++i) {
            fprintf(ofile, "%i %i %i %f %f 0\n", 
                i + 1, 
                system.atom_group[i],
                system.atom_type[i], 
                state.pos[i][0], 
                state.pos[i][1]);
        }
    }
};

class LineDumper : public Dumper {
public:
    LineDumper(const char* filename, const std::vector<std::string>& dump_names) : 
        Dumper(filename),
        dump_names(dump_names) {

    }

    void dump(const std::map<std::string, double> content, const size_t& step) {
        fprintf(ofile, "\n%zd", step);
        for (const auto& name : dump_names) {
            fprintf(ofile, "\t%f", content.at(name));
        }
    }
    
    std::vector<std::string> dump_names;
};

#endif // !ACTIVERUN_DUMPER
