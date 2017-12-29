#include "Dumper.h"

#include "System.h"
#include "Context.h"
#include <stdarg.h>

#define FORMATTER "% .6g"

void sprintf_vec(char* buffer, const Vec2& v) {
	sprintf(buffer, "% .5g % .5g", v[0], v[1]);
}

void sprintf_vec(char* buffer, const Vec3& v) {
	sprintf(buffer, "% .5g % .5g % .5g", v[0], v[1], v[2]);
}

void foprintf(FILE* file, const char* str, ...) {
	va_list args;
	va_start(args, str);
	char buffer[256];
	vsprintf(buffer, str, args);
	va_end(args);

	printf("%s", buffer);
	fprintf(file, "%s", buffer);
}

Dumper::Dumper(const char* filename, const char* mode) {
	ofile = fopen(filename, mode);
	if (!ofile) {
		fprintf(stderr, "Cannot open %s", filename);
		throw std::runtime_error("IO Error");
	}
}

void Dumper::flush() {
	fflush(ofile);
}

int Dumper::write(const char* str, ...) {

	va_list args;
	va_start(args, str);
	int ret = vfprintf(ofile, str, args);
	va_end(args);
	return ret;
}

MultiDumper::MultiDumper(const std::vector<const char*> filenames, const char* mode) {
	for (const auto& filename : filenames) {
		FILE* ofile = fopen(filename, mode);
		if (!ofile) {
			fprintf(stderr, "Cannot open %s", filename);
			throw std::runtime_error("IO Error");
		}
		ofiles.push_back(ofile);
	}
}

void MultiDumper::add_file(const char* filename, const char* mode) {
	FILE* ofile = fopen(filename, mode);
	if (!ofile) {
		fprintf(stderr, "Cannot open %s", filename);
		throw std::runtime_error("IO Error");
	}
	ofiles.push_back(ofile);
}

void MultiDumper::add_file(FILE* file) {
	ofiles.push_back(file);
}

void MultiDumper::flush_all() {
	for (const auto& ofile : ofiles) {
		fflush(ofile);
	}
}

void MultiDumper::write_all(const char* str, ...) {
	va_list args;
	va_start(args, str);
	
	for (const auto& ofile : ofiles) {
		vfprintf(ofile, str, args);
	}
	va_end(args);
}

void TrajDumper::dump(const System& system, const State& state, size_t step) {
	m_dumper.write("ITEM: TIMESTEP\n");
	m_dumper.write("%zd\n", step);
	m_dumper.write("ITEM: NUMBER OF ATOMS\n");
	m_dumper.write("%zd\n", system.atom_num);
	m_dumper.write("ITEM: BOX BOUNDS\n");
	m_dumper.write("0 %f\n", system.box[0]);
	m_dumper.write("0 %f\n", system.box[1]);
#ifdef THREE_DIMENSION
	m_dumper.write("0 %f\n", system.box[2]);
#else
	m_dumper.write("-0.25 0.25\n");
#endif // THREE_DIMENSION

	m_dumper.write("ITEM: ATOMS id mol type x y z\n");
	for (int i = 0; i < system.atom_num; ++i) {
#ifdef THREE_DIMENSION
		m_dumper.write("%i %i %i %f %f %f\n",
			i + 1,
			system.atom_group[i],
			system.atom_type[i],
			state.pos[i][0],
			state.pos[i][1],
			state.pos[i][2]);

#else
		m_dumper.write("%i %i %i %f %f 0\n",
			i + 1,
			system.atom_group[i],
			system.atom_type[i],
			state.pos[i][0],
			state.pos[i][1]);
#endif // THREE_DIMENSION

	}
}

void ThermoDumper::dump_head() {
	m_dumper.write_all("step");
	for (const auto& name : dump_names) {
		m_dumper.write_all("\t%s", name.c_str());
	}
    m_dumper.write_all("\n");
}

void ThermoDumper::dump(const std::vector<double>& value, const size_t& step) {
	m_dumper.write_all("%zd", step);
	for (const auto& v : value) {
		m_dumper.write_all("\t" FORMATTER, v);
	}
	m_dumper.write_all("\n");
}

/* Error dump */

void dump_pos(const State& state, FILE* dumpfile) {
	char buffer[256];

	fprintf(dumpfile, "Positions:\n");
	for (size_t i = 0; i < state.pos.size(); i++) {
		sprintf_vec(buffer, state.pos[i]);
		fprintf(dumpfile, "%zd    %s\n", i, buffer);
	}
}

void dump_force(const Context& context, FILE* dumpfile) {
	char buffer[256];
	fprintf(dumpfile, "Forces:\n");
	for (size_t i = 0; i < context.force_buffer[0].size(); i++) {
		fprintf(dumpfile, "%zd", i);
		for (size_t j = 0; j < context.force_buffer.size(); j++) {
			sprintf_vec(buffer, context.force_buffer[j][i]);
			fprintf(dumpfile, "       %s", buffer);
		}
		fprintf(dumpfile, "\n");
	}
}


void dump_nblist(const NeighbourList& nblist, FILE* dumpfile) {
	fprintf(dumpfile, "NeighbourList:\n");
	for (int i = 0; i < nblist.box_num[0]; i++)
		for (int j = 0; j < nblist.box_num[1]; j++) 
#ifndef THREE_DIMENSION
		{
			fprintf(dumpfile, "%d,%d: ", i, j);
			for (auto& x : nblist.at(i, j))
#else
			for(int k = 0; k < nblist.box_num[2]; k++){
				fprintf(dumpfile, "%d,%d,%d: ", i, j, k);
			for(auto& x: nblist.at(i, j, k))
#endif // !THREE_DIMENSION
			{
				fprintf(dumpfile, " %zd", x);
			}
			fprintf(dumpfile, "\n");
		}
}

void dump_snapshot(const State& state, const Context& context, const char* name) {
	FILE* dumpfile = fopen(name, "w");
	dump_pos(state, dumpfile);
	dump_force(context, dumpfile);
	dump_nblist(*context.neigh_list, dumpfile);
	fclose(dumpfile);
	fprintf(stderr, "shapshot saved to %s\n", name);
}

void dump_xyz(const State & state, const char * name)
{
    FILE* dumpfile = fopen(name, "w");
    for (auto& pos : state.pos) {
#ifndef THREE_DIMENSION
        fprintf(dumpfile, "%lf %lf\n", pos[0], pos[1]);
#else
        fprintf(dumpfile, "%lf %lf %lf\n", pos[0], pos[1], pos[2]);
#endif // !THREE_DIMENSION

    }
    fclose(dumpfile);
}

void read_xyz(State& state, const char* name) {
    FILE* dumpfile = fopen(name, "r");
    for (auto& pos : state.pos) {
#ifndef THREE_DIMENSION
        fscanf(dumpfile, "%lf%lf\n", &pos[0], &pos[1]);
#else
        fscanf(dumpfile, "%lf%lf%lf\n", &pos[0], &pos[1], &pos[2]);
#endif // !THREE_DIMENSION

    }
    fclose(dumpfile);
}