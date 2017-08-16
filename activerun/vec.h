#pragma once

#ifndef SIMUTIL_VEC_H
#define SIMUTIL_VEC_H

class Vec2 {
public:
    Vec2() {
		memset(data, 0, sizeof(double) * 2);
    }
    Vec2(double x, double y) {
        data[0] = x;
        data[1] = y;
    }

    double operator[](int i)const {
        return data[i];
    }
    double& operator[](int i) {
        return data[i];
    }

    bool operator==(const Vec2& rhs)const {
        return (data[0] == rhs.data[0]) && (data[1] == rhs.data[1]);
    }
    bool operator!= (const Vec2& rhs)const {
        return (data[0] != rhs.data[0]) || (data[1] != rhs.data[1]);
    }

    Vec2 operator+(const Vec2& rhs)const {
        return Vec2(data[0] + rhs.data[0], data[1] + rhs.data[1]);
    }
    Vec2& operator+=(const Vec2& rhs) {
        data[0] += rhs.data[0];
        data[1] += rhs.data[1];
        return *this;
    }

    Vec2 operator-()const {
        return Vec2(-data[0], -data[1]);
    }
    Vec2 operator-(const Vec2& rhs)const {
        return Vec2(data[0] - rhs.data[0], data[1] - rhs.data[1]);
    }
    Vec2& operator-=(const Vec2& rhs) {
        data[0] -= rhs.data[0];
        data[1] -= rhs.data[1];
        return *this;
    }

    // element-wise product
    Vec2 operator*(const Vec2& rhs)const {
        return Vec2(data[0] * rhs.data[0], data[1] * rhs.data[1]);
    }
    Vec2 operator*(double rhs)const {
        return Vec2(data[0] * rhs, data[1] * rhs);
    }
    Vec2& operator*=(double rhs) {
        data[0] *= rhs;
        data[1] *= rhs;
        return *this;
    }

    // element-wise division
    Vec2 operator/(const Vec2& rhs)const {
        return Vec2(data[0] / rhs.data[0], data[1] / rhs.data[1]);
    }
    Vec2 operator/(double rhs)const {
        return Vec2(data[0] / rhs, data[1] / rhs);
    }
    Vec2& operator/=(double rhs) {
        data[0] /= rhs;
        data[1] /= rhs;
        return *this;
    }
    double dot(const Vec2& rhs)const {
        return data[0] * rhs.data[0] + data[1] * rhs.data[1];
    }
    double norm2()const {
        return data[0] * data[0] + data[1] * data[1];
    }

private:
    double data[2];
};

class Vec3 {
public:
	Vec3() {
		memset(data, 0, sizeof(double) * 3);
	}
	Vec3(double x, double y, double z) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}

	double operator[](int i)const {
		return data[i];
	}
	double& operator[](int i) {
		return data[i];
	}

	bool operator==(const Vec3& rhs)const {
		return (data[0] == rhs.data[0]) && (data[1] == rhs.data[1]) && (data[2] == rhs.data[2]);
	}
	bool operator!= (const Vec3& rhs)const {
		return (data[0] != rhs.data[0]) || (data[1] != rhs.data[1]) || (data[2] != rhs.data[2]);
	}

	Vec3 operator+(const Vec3& rhs)const {
		return Vec3(data[0] + rhs.data[0], data[1] + rhs.data[1], data[2] + rhs.data[2]);
	}
	Vec3& operator+=(const Vec3& rhs) {
		data[0] += rhs.data[0];
		data[1] += rhs.data[1];
		data[2] += rhs.data[2];
		return *this;
	}

	Vec3 operator-()const {
		return Vec3(-data[0], -data[1], -data[2]);
	}
	Vec3 operator-(const Vec3& rhs)const {
		return Vec3(data[0] - rhs.data[0], data[1] - rhs.data[1], data[2] - rhs.data[2]);
	}
	Vec3& operator-=(const Vec3& rhs) {
		data[0] -= rhs.data[0];
		data[1] -= rhs.data[1];
		data[2] -= rhs.data[2];
		return *this;
	}

	// element-wise product
	Vec3 operator*(const Vec3& rhs)const {
		return Vec3(data[0] * rhs.data[0], data[1] * rhs.data[1], data[2] * rhs.data[2]);
	}
	Vec3 operator*(double rhs)const {
		return Vec3(data[0] * rhs, data[1] * rhs, data[2] * rhs);
	}
	Vec3& operator*=(double rhs) {
		data[0] *= rhs;
		data[1] *= rhs;
		data[2] *= rhs;
		return *this;
	}

	// element-wise division
	Vec3 operator/(const Vec3& rhs)const {
		return Vec3(data[0] / rhs.data[0], data[1] / rhs.data[1], data[2] / rhs.data[2]);
	}
	Vec3 operator/(double rhs)const {
		return Vec3(data[0] / rhs, data[1] / rhs, data[2] / rhs);
	}
	Vec3& operator/=(double rhs) {
		data[0] /= rhs;
		data[1] /= rhs;
		data[2] /= rhs;
		return *this;
	}
	double dot(const Vec3& rhs)const {
		return data[0] * rhs.data[0] + data[1] * rhs.data[1] + data[2] * rhs.data[2];
	}
	double norm2()const {
		return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
	}
	Vec3 cross(const Vec3& rhs)const {
		return Vec3(
			data[1] * rhs.data[2] - data[2] * rhs.data[1],
			data[2] * rhs.data[0] - data[0] * rhs.data[2],
			data[0] * rhs.data[1] - data[1] * rhs.data[0]
		);
	}

private:
	double data[2];

};

class Vec2d {
public:
    Vec2d() {
		memset(data, 0, sizeof(int) * 2);
    }
    Vec2d(const Vec2& v) {
        data[0] = (int)v[0];
        data[1] = (int)v[1];
    }
	Vec2d(int x, int y) {
		data[0] = x;
		data[1] = y;
	}
    int operator[](int i)const {
        return data[i];
    }
    int& operator[](int i) {
        return data[i];
    }
    Vec2 to_float()const {
        return Vec2((double)data[0], (double)data[1]);
    }
private:
    int data[3];
};

class Vec3d {
public:
	Vec3d() {
		memset(data, 0, sizeof(int) * 3);
	}
	Vec3d(const Vec3& v) {
		data[0] = (int)v[0];
		data[1] = (int)v[1];
		data[2] = (int)v[2];
	}
	Vec3d(int x, int y, int z) {
		data[0] = x;
		data[1] = y;
		data[2] = z;
	}
	int operator[](int i)const {
		return data[i];
	}
	int& operator[](int i) {
		return data[i];
	}
	Vec3 to_float()const {
		return Vec3((double)data[0], (double)data[1], (double)data[2]);
	}
private:
	int data[3];
};


#endif // !SIMUTIL_VEC_H

