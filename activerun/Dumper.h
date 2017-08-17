#pragma once

#ifndef ACTIVERUN_DUMPER
#define ACTIVERUN_DUMPER

#include "includes.h"
#include "System.h"
#include "Context.h"

class Dumper {
public:
    FILE* ofile;

    Dumper(const char* filename, const char* mode="w") {
        ofile = fopen(filename, mode);
        if (!ofile) {
            fprintf(stderr, "Cannot open %s", filename);
        }
    }

	Dumper(FILE* file) : ofile(file) {

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

    TrajDumper(const char* filename, bool is_restart=false) : Dumper(filename, is_restart ? "a":"w") {

    }

	void dump(const System& system, const State& state, size_t step);
};

class LineDumper : public Dumper {
public:


    LineDumper(const char* filename, const std::vector<std::string>& dump_names, bool with_output=false, bool is_restart=false) : 
        Dumper(filename, is_restart ? "a":"w"),
        dump_names(dump_names),
		with_output(with_output)
		{

    }

	void dump_head();

	void dump(const std::vector<double>& value, const size_t& step);
    
    std::vector<std::string> dump_names;
	bool with_output;
};

// used in error handling
void dump_snapshot(const State& state, const Context& context, const char* name="snapshot.txt");

void dump_xyz(const State& state, const char* name);

#endif // !ACTIVERUN_DUMPER
