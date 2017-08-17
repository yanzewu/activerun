#pragma once
#include <string>
#include <fstream>
#include "external/json.hpp"
#include "dict.h"


using json = nlohmann::json;

class Serializer {
public:

    Serializer() {
        
    }

    Serializer(const json& js) : root(js) {

    }

    void read_file(const char* filename) {

        std::ifstream infile(filename);
        if (!infile.is_open()) {
            infile.close();
            if (!root.empty()) {
                std::ofstream outfile(filename);
                outfile << root;
                outfile.close();
            }
        }
        if (root.empty()) {
            infile >> root;
        }
        else {
            json tmp;
            infile >> tmp;
            copy(tmp, root);
        }
        infile.close();
    }

    const json& at(const std::string& name)const {
        return root.at(name);
    }

    json& operator[](const std::string& name) {
        return root[name];
    }

    static void copy(const json& src, json& dst) {
        for (auto src_iter = src.begin(); src_iter != src.end(); src_iter++) {
            if (src_iter->is_object()) {
                copy(*src_iter, dst[src_iter.key()]);
            }
            else {
                dst[src_iter.key()] = *src_iter;
            }
        }
    }

    const Serializer child(const std::string& name)const {
        return Serializer(root.at(name));
    }

    template<class T>
    const T& get(const std::string& name, const T& d) {
        auto result = root.find(name);
        if (result != root.end())return result->get<T>();
        else return d;
    }
    
    Dict get_dict(const std::string& name) {
        json entry = root.at(name);
        Dict result;
        for (auto iter = entry.begin(); iter != entry.end(); iter++) {
            if (iter->is_boolean()) {
                bool val = iter.value();
                result[iter.key()] = val ? 1.0 : 0.0;
            }
            else {
                result[iter.key()] = (double)(iter.value());
            }
        }
        return result;
    }

private:

    json root;
};