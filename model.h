#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"
#include "tgaimage.h"

class Model {
private:
	std::vector<Vec3f> verts_;
	std::vector<Vec3f> tverts_;
	std::vector<Vec3f> nverts_;
	std::vector<std::vector<int> > faces_;
public:
	Model(std::string const& filename);
	~Model();
	int nverts();
	int nfaces();
	Vec3f vert(int i);
	Vec3f vert(int f, int v);
	Vec3f tvert(int i);
	Vec3f tvert(int f, int v);
	Vec3f nvert(int i);
	std::vector<int> face(int idx);
	// My addition
	TGAImage textures;
	TGAImage normals;
};

#endif //__MODEL_H__
