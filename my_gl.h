#ifndef __MY_GL_H__
#define __MY_GL_H__

#include <limits>
#include <algorithm>
#include "geometry.h"
#include "tgaimage.h"


extern Matrix ViewPort;
extern Matrix ModelView;
extern Matrix Projection;

TGAColor get_texture(Vec3f *tv, Vec3f bary, TGAImage &texture);
Vec3f get_normal(Vec3f *nv, Vec3f bary, TGAImage &normalmap);

Vec3f barycentric(Vec2f *pts, Vec2f P);
float max_elevation_angle(float *z, Vec2f pt, Vec2f dir, TGAImage &image);
void triangle_blank(Vec3f *v,float *zbuffer, TGAImage& image);
Matrix view_frame(Vec3f eye, Vec3f centre, Vec3f up);
Matrix perspective(float c);
Matrix viewport(int x, int y, int w, int h, int depth);
void fill_shadow_buffer(Vec3f *v, float *sb, TGAImage & image);

struct Shader{
	virtual ~Shader() {};
	virtual Vec3f vertex(int face, int vert) = 0;
	virtual void fragment(Vec3f bary, TGAColor &c, Matrix &MNinv,TGAImage &image) = 0;
};


void triangle(Vec3f v[3], Shader &shader, Matrix &M, float *zbuffer, float *sbuffer, TGAImage &image);

#endif
