#pragma once

#ifndef SIMUTIL_NEIGHLIST_H
#define SIMUTIL_NEIGHlIST_H


#include "vec.h"
#include <vector>

#ifdef THREE_DIMENSION
#define NeighbourList NeighbourList3
#else
#define NeighbourList NeighbourList2
#endif // THREE_DIMENSION



class NeighbourList2 {
public:

    NeighbourList2(const Vec2& box, double cutoff, bool using_ghost=false) {
		printf("\nInitializing celllist\n\n");

        box_num = Vec2d(box / cutoff);
		printf("Total %d x %d boxes\n", box_num[0], box_num[1]);

        unit = box / box_num.to_float();
		printf("Box size=%.4f, %.4f\n", unit[0], unit[1]);

		this->using_ghost = using_ghost;

        //2d
		if (!using_ghost) {
			cells.resize(box_num[0] * box_num[1]);
			box_num_real = box_num;
		}
		else {
			printf("Using ghost\n");
			cells.resize((box_num[0] + 2) * (box_num[1] + 2));
			box_num_real = Vec2d(box_num[0] + 2, box_num[1] + 2);
		}
    }

    void build_from_pos(const std::vector<Vec2>& pos) {
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

    Vec2 unit;
    Vec2d box_num;
	Vec2d box_num_real;
	bool using_ghost;

    std::vector<std::vector<size_t> > cells;
};

class NeighbourList3 {
public:

	NeighbourList3(const Vec3& box, double cutoff, bool using_ghost = false) {
		printf("\nInitializing celllist\n\n");

		box_num = Vec3d(box / cutoff);
		printf("Total %d x %d x %d boxes\n", box_num[0], box_num[1], box_num[2]);

		unit = box / box_num.to_float();
		printf("Box size=%.4f, %.4f, %.4f\n", unit[0], unit[1], unit[2]);

		this->using_ghost = using_ghost;

		//2d
		if (!using_ghost) {
			cells.resize(box_num[0] * box_num[1]);
			box_num_real = box_num;
		}
		else {
			printf("Using ghost\n");
			cells.resize((box_num[0] + 2) * (box_num[1] + 2) * (box_num[2] + 2));
			box_num_real = Vec3d(box_num[0] + 2, box_num[1] + 2, box_num[2] + 2);
		}
	}

	void build_from_pos(const std::vector<Vec3>& pos) {
		for (auto& cell : cells) {
			cell.clear();
		}
		if (using_ghost) {
			for (size_t i = 0; i < pos.size(); i++) {
				at((int)(pos[i][0] / unit[0]) + 1, 
					(int)(pos[i][1] / unit[1]) + 1, 
					(int)(pos[i][2] / unit[2]) + 1).push_back(i);
			}
			synchronize_ghost();
		}
		else {
			for (size_t i = 0; i < pos.size(); i++) {
				at((int)(pos[i][0] / unit[0]), 
					(int)(pos[i][1] / unit[1]),
					(int)(pos[i][1] / unit[2])).push_back(i);
			}
		}
	}

	void synchronize_ghost() {

		// xy - plane
		for (int i = 1; i <= box_num[0]; i++)
		for (int j = 1; j <= box_num[1]; j++){
			at(i, j, 0) = at(i, j, box_num[2]);
			at(i, j, box_num[2] + 1) = at(i, j, 1);
		}

		// yz - plane
		for (int i = 1; i <= box_num[1]; i++)
		for (int j = 1; j <= box_num[2]; j++){
			at(0, i, j) = at(box_num[0], i, j);
			at(box_num[0] + 1, i, j) = at(1, i, j);
		}

		// zx - plane
		for (int i = 1; i <= box_num[2]; i++)
		for (int j = 1; j <= box_num[0]; j++) {
			at(j, 0, i) = at(j, box_num[1], i);
			at(j, box_num[1] + 1, i) = at(j, 1, i);
		}

		at(0, 0, 0) = at(box_num[0], box_num[1], box_num[2]);
		at(box_num[0] + 1, 0,              0             ) = at(1,          box_num[1], box_num[2]);
		at(0,              box_num[1] + 1, 0             ) = at(box_num[0], 1,          box_num[2]);
		at(0,              0,              box_num[2] + 1) = at(box_num[0], box_num[1], 1);
		at(box_num[0] + 1, box_num[1] + 1, 0             ) = at(1,          1,          box_num[2]);
		at(box_num[0] + 1, 0,              box_num[2] + 1) = at(1,          box_num[1], 1);
		at(0,              box_num[1] + 1, box_num[2] + 1) = at(box_num[0], 1,          1);
		at(box_num[0] + 1, box_num[1] + 1, box_num[2] + 1) = at(1,          1,          1);
	}

	inline std::vector<size_t>& at(int x, int y, int z) {
		return cells[(x * box_num_real[1] + y) * box_num_real[2] + z];
	}
	// no pbc wrapper
	inline const std::vector<size_t>& at(int x, int y, int z)const {
		return cells[(x * box_num_real[1] + y) * box_num_real[2] + z];
	}
	// with pbc wrapper
	inline const std::vector<size_t>& operator()(int x, int y, int z)const {
		x = (x >= 0 ? x % box_num_real[0] : x + box_num_real[0]);
		y = (y >= 0 ? y % box_num_real[1] : y + box_num_real[1]);
		z = (z >= 0 ? z % box_num_real[2] : z + box_num_real[2]);
		return at(x, y, z);
	}

	Vec3 unit;
	Vec3d box_num;
	Vec3d box_num_real;
	bool using_ghost;

	std::vector<std::vector<size_t> > cells;
};


#endif // !SIMUTIL_NEIGHLIST_H

