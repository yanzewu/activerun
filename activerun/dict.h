#pragma once

#include <map>
#include <string>

class Dict {
public:

	Dict() {

	}
	Dict(const std::initializer_list<std::pair<const std::string, double> >& init_list) : data(init_list) {
		
	}
	double& operator[](const std::string& s) {
		return data[s];
	}
	const double& operator[](const std::string& s)const {
		return data.at(s);
	}
	double at(const std::string& s)const {
		try {
			return data.at(s);
		}
		catch (const std::out_of_range&) {
			fprintf(stderr, "Cannot initialize essential element: %s\n", s.c_str());
			throw;
		}
	}
	double get(const std::string& s, double d)const {
		auto r = data.find(s);
		if (r == data.end())return d;
		else return r->second;
	}

protected:

	std::map<std::string, double> data;
};