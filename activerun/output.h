#pragma once

#ifndef OUTPUT_H
#define OUTPUT_H

#include <stdio.h>
#include <vector>
#include <stdarg.h>

class Dumper {
public:
	FILE* ofile;

	Dumper() : ofile(NULL) {

	}

	Dumper(const char* filename, const char* mode = "w");

	Dumper(FILE* file) : ofile(file) {

	}
	Dumper(const Dumper& d) : ofile(d.ofile) {

	}

	~Dumper() {
		if (ofile) fclose(ofile);
		ofile = NULL;
	}

	void flush();

	int write(const char* str, ...);

};

class MultiDumper {
public:

	MultiDumper(const std::vector<const char*> filenames, const char* mode = "w");

	MultiDumper(const std::vector<FILE*>& files) : ofiles(files) {

	}

	MultiDumper(const MultiDumper& md) : ofiles(md.ofiles) {

	}

	~MultiDumper() {
		for (auto& ofile : ofiles) {
			if (ofile) fclose(ofile);
			ofile = NULL;
		}
	}

	void add_file(const char* filename, const char* mode = "w");

	void add_file(FILE* file);

	void flush_all();

	void write_all(const char* str, ...);

private:
	std::vector<FILE*> ofiles;
};

#endif