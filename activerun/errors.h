#pragma once

#include "System.h"
#include "Context.h"


void dump_pos(const State& state, FILE* dumpfile) {
	fprintf(dumpfile, "Positions:\n");
	for (size_t i = 0; i < state.pos.size(); i++) {
		fprintf(dumpfile, "%zd    % .5e,%.5e\n", i, state.pos[i][0], state.pos[i][1]);
	}
}

void dump_force(const Context& context, FILE* dumpfile) {
	fprintf(dumpfile, "Forces:\n");
	for (size_t i = 0; i < context.force_buffer[0].size(); i++) {
		fprintf(dumpfile, "%zd", i);
		for (size_t j = 0; j < context.force_buffer.size(); j++) {
			fprintf(dumpfile, "    % .5e,% .5e", context.force_buffer[j][i][0], context.force_buffer[j][i][1]);
		}
		fprintf(dumpfile, "\n");
	}
}


void dump_nblist(const NeighbourList& nblist, FILE* dumpfile) {
	fprintf(dumpfile, "NeighbourList:\n");
	for(size_t i = 0; i < nblist.box_num[0]; i++)
		for (size_t j = 0; j < nblist.box_num[1]; j++) {
			fprintf(dumpfile, "%d,%d  ", i, j);
			for (auto& x : nblist.at(i, j)) {
				fprintf(dumpfile, " %d", x);
			}
			fprintf(dumpfile, "\n");
		}
}

void dump(const State& state, const Context& context) {
	FILE* dumpfile = fopen("error.log", "w");
	dump_pos(state, dumpfile);
	dump_force(context, dumpfile);
	fclose(dumpfile);
}