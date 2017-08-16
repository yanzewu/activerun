#include "Dumper.h"

#pragma once

#include "System.h"
#include "Context.h"
#include <stdarg.h>

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

void TrajDumper::dump(const System& system, const State& state, size_t step) {
	fprintf(ofile, "ITEM: TIMESTEP\n");
	fprintf(ofile, "%zd\n", step);
	fprintf(ofile, "ITEM: NUMBER OF ATOMS\n");
	fprintf(ofile, "%zd\n", system.atom_num);
	fprintf(ofile, "ITEM: BOX BOUNDS\n");
	fprintf(ofile, "0 %f\n", system.box[0]);
	fprintf(ofile, "0 %f\n", system.box[1]);
#ifdef THREE_DIMENSION
	fprintf(ofile, "0 %f\n", system.box[2]);
#else
	fprintf(ofile, "-0.25 0.25\n");
#endif // THREE_DIMENSION

	fprintf(ofile, "ITEM: ATOMS id mol type x y z\n");
	for (int i = 0; i < system.atom_num; ++i) {
#ifdef THREE_DIMENSION
		fprintf(ofile, "%i %i %i %f %f %f\n",
			i + 1,
			system.atom_group[i],
			system.atom_type[i],
			state.pos[i][0],
			state.pos[i][1],
			state.pos[i][2]);

#else
		fprintf(ofile, "%i %i %i %f %f 0\n",
			i + 1,
			system.atom_group[i],
			system.atom_type[i],
			state.pos[i][0],
			state.pos[i][1]);
#endif // THREE_DIMENSION

	}
}

void LineDumper::dump_head() {
	if (with_output) {
		foprintf(ofile, "step");
		for (const auto& name : dump_names) {
			foprintf(ofile, "\t%s", name.c_str());
		}
	}
	else {
		fprintf(ofile, "step");
		for (const auto& name : dump_names) {
			fprintf(ofile, "\t%s", name.c_str());
		}
	}
}

void LineDumper::dump(const std::vector<double>& value, const size_t& step) {
	if (with_output) {
		foprintf(ofile, "\n%zd", step);
		for (const auto& v : value) {
			foprintf(ofile, "\t%f", v);
		}
	}
	else {
		fprintf(ofile, "\n%zd", step);
		for (const auto& v : value) {
			fprintf(ofile, "\t%f", v);
		}

	}
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
		for (int j = 0; j < nblist.box_num[1]; j++) {
			fprintf(dumpfile, "%d,%d  ", i, j);
			for (auto& x : nblist.at(i, j)) {
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