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

template<class T>
inline void copy_vec(std::vector<T>& dst, const std::vector<T>& src) {
    dst.assign(src.begin(), src.end());
}

class NeighbourList2 {
public:

    NeighbourList2(const Vec2& box, double cutoff, bool using_ghost=false) {
        box_num = Vec2d(box / cutoff);
        unit = box / box_num.to_float();

		this->using_ghost = using_ghost;

        //2d
		if (!using_ghost) {
			cells.resize(box_num[0] * box_num[1]);
			box_num_real = box_num;
		}
		else {
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
		copy_vec(at(0, 0), at(box_num[0], box_num[1]));
        copy_vec(at(0, box_num[1] + 1), at(box_num[0], 1));
        copy_vec(at(box_num[0] + 1, 0), at(1, box_num[1]));
        copy_vec(at(box_num[0] + 1, box_num[1] + 1), at(1, 1));
	}

    void compute_ghost_pos(std::vector<size_t>& ghost_id, std::vector<Vec2d>& ghost_offset)const {

        auto append_ghost = [&](const std::vector<size_t>& box, const Vec2d& offset, std::vector<size_t>& gid, std::vector<Vec2d>& goff) {
            gid.insert(gid.end(), box.begin(), box.end());
            goff.insert(goff.end(), box.size(), offset);
        };

        for (int i = 1; i <= box_num[0]; i++) {
            append_ghost(at(i, 0), { 0, -1 }, ghost_id, ghost_offset);
            append_ghost(at(i, box_num[1] + 1), { 0, 1 }, ghost_id, ghost_offset);
        }
        for (int i = 1; i <= box_num[1]; i++) {
            append_ghost(at(0, i), { -1, 0 }, ghost_id, ghost_offset);
            append_ghost(at(box_num[0] + 1, i), { 1, 0 }, ghost_id, ghost_offset);
        }

        append_ghost(at(0, 0), { -1, -1 }, ghost_id, ghost_offset);
        append_ghost(at(0, box_num[1] + 1), { -1, 1 }, ghost_id, ghost_offset);
        append_ghost(at(box_num[0] + 1, 0), { 1, -1 }, ghost_id, ghost_offset);
        append_ghost(at(box_num[0] + 1, box_num[1] + 1), { 1, 1 }, ghost_id, ghost_offset);
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
		box_num = Vec3d(box / cutoff);
        if (using_ghost && (box_num[0] == 1 || box_num[1] == 1 || box_num[2] == 1)) {
            fprintf(stderr, "Error: box too small\n");
            throw std::runtime_error("Cannot initialize neighbourlist");
        }
		unit = box / box_num.to_float();

		this->using_ghost = using_ghost;

		//2d
		if (!using_ghost) {
			cells.resize(box_num[0] * box_num[1] * box_num[2]);
			box_num_real = box_num;
		}
		else {
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
            copy_vec(at(i, j, 0), at(i, j, box_num[2]));
            copy_vec(at(i, j, box_num[2] + 1), at(i, j, 1));
		}

		// yz - plane
		for (int i = 1; i <= box_num[1]; i++)
		for (int j = 1; j <= box_num[2]; j++){
            copy_vec(at(0, i, j), at(box_num[0], i, j));
            copy_vec(at(box_num[0] + 1, i, j), at(1, i, j));
		}

		// zx - plane
		for (int i = 1; i <= box_num[2]; i++)
		for (int j = 1; j <= box_num[0]; j++) {
            copy_vec(at(j, 0, i), at(j, box_num[1], i));
            copy_vec(at(j, box_num[1] + 1, i), at(j, 1, i));
		}

		// z axis
		for (int i = 1; i <= box_num[2]; i++) {
            copy_vec(at(0,              0,              i), at(box_num[0], box_num[1], i));
            copy_vec(at(0,              box_num[1] + 1, i), at(box_num[0], 1,          i));
            copy_vec(at(box_num[0] + 1, 0,              i), at(1,          box_num[1], i));
            copy_vec(at(box_num[0] + 1, box_num[1] + 1, i), at(1,          1,          i));
		}

		// x axis
		for (int i = 1; i <= box_num[0]; i++) {
            copy_vec(at(i, 0,              0             ), at(i, box_num[1], box_num[2]));
            copy_vec(at(i, 0,              box_num[2] + 1), at(i, box_num[1], 1));
            copy_vec(at(i, box_num[1] + 1, 0             ), at(i, 1,          box_num[2]));
            copy_vec(at(i, box_num[1] + 1, box_num[2] + 1), at(i, 1,          1));
		}

		// y axis
		for (int i = 1; i <= box_num[1]; i++) {
            copy_vec(at(0,              i, 0             ), at(box_num[0], i, box_num[2]));
            copy_vec(at(box_num[0] + 1, i, 0             ), at(1,          i, box_num[2]));
            copy_vec(at(0             , i, box_num[2] + 1), at(box_num[0], i, 1));
            copy_vec(at(box_num[0] + 1, i, box_num[2] + 1), at(1,          i, 1));
		}

        copy_vec(at(0, 0, 0), at(box_num[0], box_num[1], box_num[2]));
        copy_vec(at(box_num[0] + 1, 0,              0             ), at(1,          box_num[1], box_num[2]));
        copy_vec(at(0,              box_num[1] + 1, 0             ), at(box_num[0], 1,          box_num[2]));
        copy_vec(at(0,              0,              box_num[2] + 1), at(box_num[0], box_num[1], 1));
        copy_vec(at(box_num[0] + 1, box_num[1] + 1, 0             ), at(1,          1,          box_num[2]));
        copy_vec(at(box_num[0] + 1, 0,              box_num[2] + 1), at(1,          box_num[1], 1));
        copy_vec(at(0,              box_num[1] + 1, box_num[2] + 1), at(box_num[0], 1,          1));
        copy_vec(at(box_num[0] + 1, box_num[1] + 1, box_num[2] + 1), at(1,          1,          1));
	}

    void compute_ghost_pos(std::vector<size_t>& ghost_id, std::vector<Vec3d>& ghost_offset)const {

        auto append_ghost = [&](const std::vector<size_t>& box, const Vec3d& offset, std::vector<size_t>& gid, std::vector<Vec3d>& goff) {
            gid.insert(gid.end(), box.begin(), box.end());
            goff.insert(goff.end(), box.size(), offset);
        };
        // xy - plane
        for (int i = 1; i <= box_num[0]; i++)
        for (int j = 1; j <= box_num[1]; j++) {
            append_ghost(at(i, j, 0), { 0, 0, -1 }, ghost_id, ghost_offset);
            append_ghost(at(i, j, box_num[2] + 1), { 0, 0, 1 }, ghost_id, ghost_offset);
        }

        // yz - plane
        for (int i = 1; i <= box_num[1]; i++)
        for (int j = 1; j <= box_num[2]; j++) {
            append_ghost(at(0, i, j), { -1, 0, 0 }, ghost_id, ghost_offset);
            append_ghost(at(box_num[0] + 1, i, j), { 1, 0, 0 }, ghost_id, ghost_offset);
        }

        // zx - plane
        for (int i = 1; i <= box_num[2]; i++)
        for (int j = 1; j <= box_num[0]; j++) {
            append_ghost(at(j, 0, i), { 0, -1, 0 }, ghost_id, ghost_offset);
            append_ghost(at(j, box_num[1] + 1, i), { 0, 1, 0 }, ghost_id, ghost_offset);
        }

        // z axis
        for (int i = 1; i <= box_num[2]; i++) {
            append_ghost(at(0, 0, i), { -1, -1, 0 }, ghost_id, ghost_offset);
            append_ghost(at(0, box_num[1] + 1, i), { -1, 1, 0 }, ghost_id, ghost_offset);
            append_ghost(at(box_num[0] + 1, 0, i), { 1, -1, 0 }, ghost_id, ghost_offset);
            append_ghost(at(box_num[0] + 1, box_num[1] + 1, i), { 1, 1, 0 }, ghost_id, ghost_offset);
        }

        // x axis
        for (int i = 1; i <= box_num[0]; i++) {
            append_ghost(at(i, 0, 0), { 0, -1, -1 }, ghost_id, ghost_offset);
            append_ghost(at(i, 0, box_num[2] + 1), { 0, -1, 1 }, ghost_id, ghost_offset);
            append_ghost(at(i, box_num[1] + 1, 0), { 0, 1, -1 }, ghost_id, ghost_offset);
            append_ghost(at(i, box_num[1] + 1, box_num[2] + 1), { 0, 1, 1 }, ghost_id, ghost_offset);
        }

        // y axis
        for (int i = 1; i <= box_num[1]; i++) {
            append_ghost(at(0, i, 0), { -1, 0, -1 }, ghost_id, ghost_offset);
            append_ghost(at(box_num[0] + 1, i, 0), { 1, 0, -1 }, ghost_id, ghost_offset);
            append_ghost(at(0, i, box_num[2] + 1), { -1, 0, 1 }, ghost_id, ghost_offset);
            append_ghost(at(box_num[0] + 1, i, box_num[2] + 1), { 1, 0, 1 }, ghost_id, ghost_offset);
        }

        append_ghost(at(0, 0, 0), { -1, -1, -1 }, ghost_id, ghost_offset);
        append_ghost(at(box_num[0] + 1, 0, 0), { 1, -1, -1 }, ghost_id, ghost_offset);
        append_ghost(at(0, box_num[1] + 1, 0), { -1, 1, -1 }, ghost_id, ghost_offset);
        append_ghost(at(0, 0, box_num[2] + 1), { -1, -1, 1 }, ghost_id, ghost_offset);
        append_ghost(at(box_num[0] + 1, box_num[1] + 1, 0), { 1, 1, -1 }, ghost_id, ghost_offset);
        append_ghost(at(box_num[0] + 1, 0, box_num[2] + 1), { 1, -1, 1 }, ghost_id, ghost_offset);
        append_ghost(at(0, box_num[1] + 1, box_num[2] + 1), { -1, 1, 1 }, ghost_id, ghost_offset);
        append_ghost(at(box_num[0] + 1, box_num[1] + 1, box_num[2] + 1), { 1, 1, 1 }, ghost_id, ghost_offset);

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

