#pragma once

#ifndef ACTIVERUN_SYSTEM_H
#define ACTIVERUN_SYSTEM_H
#include <map>

#include "includes.h"
#include "Input.h"



/* State for i/o, including:
pos, stress, energy, temperature, etc.
*/
class State {
public:

    std::vector<Vec> pos;

    void init(const DataFile& datafile) {
#ifndef THREE_DIMENSION
		pos.reserve(datafile.pos.size());
		for (const auto& p : datafile.pos) {
			pos.emplace_back(Vec2(p[0], p[1]));
		}
#else
        pos = datafile.pos;
#endif // !THREE_DIMENSION
    }

	void write_data(DataFile& datafile) {
#ifndef THREE_DIMENSION
		datafile.pos.resize(0);
		for (const auto& p : pos) {
			datafile.pos.emplace_back(Vec3(p[0], p[1], 0));
		}
#else
		datafile.pos = pos;

#endif // !THREE_DIMENSION
	}
};


/* stores everything of system other than kinetics (positions, forces, velocities)

*/
class System {
public:

    void init(const DataFile& datafile) {
        atom_num = datafile.atom_num;
        atom_group = datafile.atom_group;
        atom_type = datafile.atom_type;
        atom_mass.resize(atom_num);
        atom_attributes.resize(datafile.type_pair_coeff[0].size());
        for (auto& attr : atom_attributes) {
            attr.resize(atom_num);
        }

        for (size_t i = 0; i < atom_num; i++) {
            int mytype = atom_type[i] - 1;

            atom_mass[i] = datafile.type_mass[mytype];
            for (size_t k = 0; k < atom_attributes.size(); k++) {
                atom_attributes[k][i] = datafile.type_pair_coeff[mytype][k];
            }
        }
#ifndef THREE_DIMENSION

		box = Vec2(datafile.box[0], datafile.box[1]);
#else
		box = datafile.box;
#endif // !THREE_DIMENSION
    }

    const std::vector<double>& get_attr(const std::string& name)const {
        return atom_attributes.at(attribute_names.at(name));
    }
    std::vector<double>& set_attr(const std::string& name) {
        return atom_attributes[attribute_names.at(name)];
    }
    std::vector<double>& add_attr(const std::string& name) {
        atom_attributes.push_back(std::vector<double>());
        atom_attributes.back().resize(atom_num);
        attribute_names[name] = (int)atom_attributes.size() - 1;
        return atom_attributes.back();
    }
    void set_name(const std::string& name, int index) {
        attribute_names[name] = index;
    }


    size_t atom_num;

    std::vector<int> atom_group;    // group of each item
    std::vector<int> atom_type;     // type of each item
    std::vector<double> atom_mass;
	std::vector<std::vector<double> > atom_attributes;
	std::map<std::string, int> attribute_names;


    Vec2 box;

};


#endif // !ACTIVERUN_SYSTEM_H

