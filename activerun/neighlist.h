#pragma once

#ifndef SIMUTIL_NEIGHLIST_H
#define SIMUTIL_NEIGHlIST_H


#include "vec.h"
#include <vector>

#define Vec Vec2
#define Vecd Vec2d

class NeighbourList {
public:

    NeighbourList(const Vec& box, double cutoff, bool using_ghost=false) {
        box_num = Vecd(box / cutoff);

        unit = box / box_num.to_float();

		this->using_ghost = using_ghost;

        //2d
		if (!using_ghost) {
			cells.resize(box_num[0] * box_num[1]);
			box_num_real = box_num;
		}
		else {
			cells.resize((box_num[0] + 2) * (box_num[1] + 2));
			box_num_real = Vecd(box_num[0] + 2, box_num[1] + 2);
		}
    }

    void build_from_pos(const std::vector<Vec>& pos) {
        for (auto& cell : cells) {
            cell.clear();
        }
		if (using_ghost) {
			for (size_t i = 0; i < pos.size(); i++) {
				at((int)(pos[i][0] / unit[0]) + 1, (int)(pos[i][1] / unit[1]) + 1).push_back(i);
			}
			synchronize_ghost();
		}
		else {
			for (size_t i = 0; i < pos.size(); i++) {
				at((int)(pos[i][0] / unit[0]), (int)(pos[i][1] / unit[1])).push_back(i);
			}
		}
    }

	void synchronize_ghost() {
		for (int i = 1; i <= box_num[0]; i++) {
			at(i, 0) = at(i, box_num[1]);
			at(i, box_num[1] + 1) = at(i, 1);
		}
		for (int i = 1; i <= box_num[1]; i++) {
			at(0, i) = at(box_num[0], i);
			at(box_num[0] + 1, i) = at(1, i);
		}
		at(0, 0) = at(box_num[0], box_num[1]);
		at(0, box_num[1] + 1) = at(box_num[0], 1);
		at(box_num[0] + 1, 0) = at(1, box_num[1]);
		at(box_num[0] + 1, box_num[1] + 1) = at(1, 1);
	}

    inline std::vector<size_t>& at(int x, int y) {
//		if (y >= box_num_real[1])throw std::out_of_range(std::to_string(y));
        return cells[x * box_num_real[1] + y];
    }
    // no pbc wrapper
    inline const std::vector<size_t>& at(int x, int y)const {
        // 2d
        return cells[x * box_num_real[1] + y];
    }
    // with pbc wrapper
    inline const std::vector<size_t>& operator()(int x, int y)const {
        // 2d
        x = (x >= 0 ? x % box_num_real[0] : x + box_num_real[0]);
        y = (y >= 0 ? y % box_num_real[1] : y + box_num_real[1]);

        return at(x, y);
    }

    Vec unit;
    Vecd box_num;
	Vecd box_num_real;
	bool using_ghost;

    std::vector<std::vector<size_t> > cells;
};


#endif // !SIMUTIL_NEIGHLIST_H

