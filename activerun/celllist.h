#pragma once

#include "arrayutil.h"


#define Vec Vec2
#define Vecd Vec2d

class NeighbourList {
public:

    NeighbourList(const Vec& box, double cutoff) {
        box_num = Vecd(box / cutoff);

        unit = box / box_num.to_float();

        //2d
        cells.resize(box_num[0] * box_num[1]);
    }

    void build_from_pos(const std::vector<Vec>& pos) {
        for (auto& cell : cells) {
            cell.clear();
        }
        for (size_t i = 0; i < pos.size(); i++) {
            at(pos[i][0] / unit[0], pos[i][1] / unit[1]).push_back(i);
        }
    }
    inline std::vector<size_t>& at(int x, int y) {
        return cells[x * box_num[1] + y];
    }
    // no pbc wrapper
    inline const std::vector<size_t>& at(int x, int y)const {
        // 2d
        return cells[x * box_num[1] + y];
    }
    // with pbc wrapper
    inline const std::vector<size_t>& operator()(int x, int y)const {
        // 2d
        x = (x >= 0 ? x % box_num[0] : x + box_num[0]);
        y = (y >= 0 ? y % box_num[1] : y + box_num[1]);

        return at(x, y);
    }

    Vec unit;
    Vecd box_num;

    std::vector<std::vector<size_t> > cells;
};
