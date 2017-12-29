#pragma once

#ifndef ACTIVERUN_DUMPER
#define ACTIVERUN_DUMPER

#include "includes.h"
#include "output.h"
#include "System.h"
#include "Context.h"

class TrajDumper {
public:

	TrajDumper(const char* filename, bool is_restart = false) : m_dumper(filename, is_restart ? "a" : "w") {

	}

	void dump(const System& system, const State& state, size_t step);

	void flush();
private:
	Dumper m_dumper;
};


class ThermoDumper {
public:


	ThermoDumper(const char* filename, const MultiDumper& logfile, const std::vector<std::string>& dump_names, bool is_restart=false) :
		m_dumper(logfile),
        dump_names(dump_names)
		{
		m_dumper.add_file(filename, is_restart ? "a" : "w");
    }

	void dump_head();

	void dump(const std::vector<double>& value, const size_t& step);
    
    std::vector<std::string> dump_names;
	
private:
	MultiDumper m_dumper;
};

// used in error handling
void dump_snapshot(const State& state, const Context& context, const char* name="snapshot.txt");

void dump_xyz(const State& state, const char* name);

#endif // !ACTIVERUN_DUMPER
