#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include <cmath>
#include <vector>

class Matrix;
template <class t> struct Vec2;
template <class t> struct Vec3;
template <class t> struct Vec4;
typedef Vec2<float> Vec2f;
typedef Vec2<int>   Vec2i;
typedef Vec3<float> Vec3f;
typedef Vec4<float> Vec4f;
typedef Vec3<int>   Vec3i;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class t> struct Vec2 {
	union {
		struct {t u, v;};
		struct {t x, y;};
		t raw[2];
	};
	Vec2() : u(0), v(0) {}
	Vec2(t _u, t _v) : u(_u),v(_v) {}
	inline Vec2<t> operator +(const Vec2<t> &V) const { return Vec2<t>(u+V.u, v+V.v); }
	inline Vec2<t> operator -(const Vec2<t> &V) const { return Vec2<t>(u-V.u, v-V.v); }
	inline Vec2<t> operator *(float f)          const { return Vec2<t>(u*f, v*f); }
	float norm () const{ return std::sqrt(u*u+v*v); }
	inline t operator[](const int i) { return raw[i]; }
	template <class > friend std::ostream& operator<<(std::ostream& s, Vec2<t>& v);
};

template <class t> struct Vec3 {
	union {
		struct {t x, y, z;};
		struct { t ivert, iuv, inorm; };
		t raw[3];
	};
	Vec3() : x(0), y(0), z(0) {}
	Vec3(t _x, t _y, t _z) : x(_x),y(_y),z(_z) {}
	Vec3(const Matrix m);
	template <class u> Vec3<t>(const Vec3<u> &v);
	template <class u> Vec3<t>(const Vec4<u> &v);
	inline Vec3<t> operator +(const Vec3<t> &v) const { return Vec3<t>(x+v.x, y+v.y, z+v.z); }
	inline Vec3<t> operator +=(const Vec3<t> &v) const { x+=v.x;y+=v.y;z+=v.z;return *this; }
	inline Vec3<t> operator -(const Vec3<t> &v) const { return Vec3<t>(x-v.x, y-v.y, z-v.z); }
	inline Vec3<t> operator *(float f)          const { return Vec3<t>(x*f, y*f, z*f); }
	inline Vec3<t> operator ^(const Vec3<t> &v) const {

	   	return Vec3<t>(y*v.z-z*v.y, z*v.x-x*v.z, x*v.y-y*v.x); 
	}
	inline t       operator *(const Vec3<t> &v) const {
	   	return x*v.x + y*v.y + z*v.z; 
	}
	float norm () const { return std::sqrt(x*x+y*y+z*z); }
	inline t operator[](const int i) { return raw[i]; }
	Vec3<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
	template <class > friend std::ostream& operator<<(std::ostream& s, Vec3<t>& v);
};

template <class t> struct Vec4 {
	union {
		struct {t x, y, z, c;};
		struct { t ivert, iuv, inorm, w; };
		t raw[4];
	};
	Vec4() : x(0), y(0), z(0), c(1) {}
	Vec4(t _x, t _y, t _z, t _c) : x(_x),y(_y),z(_z), c(_c) {}
	Vec4(Matrix m);
	template <class u> Vec4<t>(const Vec4<u> &v);
	template <class u> Vec4<t>(const Vec3<u> &v);
	inline Vec4<t> operator +(const Vec4<t> &v) const { return Vec4<t>(x+v.x, y+v.y, z+v.z, c+v.c); }
	inline Vec4<t> operator -(const Vec4<t> &v) const { return Vec4<t>(x-v.x, y-v.y, z-v.z, c-v.c); }
	inline Vec4<t> operator *(float f)          const { return Vec4<t>(x*f, y*f, z*f, c*f); }
	inline t       operator *(const Vec4<t> &v) const { return x*v.x + y*v.y + z*v.z + c*v.c; }
	float norm () const { return std::sqrt(x*x+y*y+z*z+c*c); }
	inline t operator[](const int i) { return raw[i]; }
	Vec3<t> & normalize(t l=1) { *this = (*this)*(l/norm()); return *this; }
	template <class > friend std::ostream& operator<<(std::ostream& s, Vec4<t>& v);
};

template <class t> std::ostream& operator<<(std::ostream& s, Vec2<t>& v) {
	s << "(" << v.x << ", " << v.y << ")\n";
	return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec3<t>& v) {
	s << "(" << v.x << ", " << v.y << ", " << v.z << ")\n";
	return s;
}

template <class t> std::ostream& operator<<(std::ostream& s, Vec4<t>& v) {
	s << "(" << v.x << ", " << v.y << ", " << v.z << ", " << v.c <<  ")\n";
	return s;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Matrix{
	std::vector<std::vector<float>> m;
	int rows, cols;
public:
	Matrix( int r = 4, int c = 4);
	Matrix(Vec3f v);
	Matrix(Vec4f v);
	int nrows();
	int ncols();
	static Matrix identity(int dimensions);
	std::vector<float>& operator[](const int i);
	Matrix operator*(const Matrix& a);
	Matrix transpose();
	Matrix inverse();
	// TODO: Make it so that this only gets computed once, not every time the shader is run.
	inline Matrix inverse_transpose() {return transpose().inverse(); }
	friend std::ostream& operator<<(std::ostream& s, Matrix& m);
};

#endif //__GEOMETRY_H__
