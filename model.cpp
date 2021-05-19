#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(std::string const& filename) : verts_(), tverts_(), nverts_(), faces_(), textures(), normals() {
    std::ifstream in;
    in.open ((filename+".obj").c_str(), std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } else if (!line.compare(0,3, "vt ")) {
            iss >> trash >> trash >> trash;
            Vec3f v;
            for (int i=0;i<3;i++){
				iss >> v.raw[i];
			}
            tverts_.push_back(v);
        } else if (!line.compare(0,3, "vn ")) {
            iss >> trash >> trash >> trash ;
            Vec3f v;
            for (int i=0;i<3;i++){
				iss >> v.raw[i];
			}
            nverts_.push_back(v);
		} else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            int idx1, idx2, idx3;
            iss >> trash;
			// while the lines parse like "int1/int2/int3" store int in idx;
            while (iss >> idx1 >> trash >> idx2 >> trash >> idx3) {
                idx1--; // in wavefront obj all indices start at 1, not zero
                idx2--;
                idx3--; 
                f.push_back(idx1);
                f.push_back(idx2);
                f.push_back(idx3); // f should be vector of length 9
            }
            faces_.push_back(f);
		}
    }
    std::cerr << "# v# " << verts_.size() << " vt# " << tverts_.size() << " f# "  << faces_.size() << std::endl;
	textures.read_tga_file(filename + "_diffuse.tga");
	normals.read_tga_file(filename + "_nm.tga");
	textures.flip_vertically();
	normals.flip_vertically();
}

Model::~Model() {
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

Vec3f Model::vert(int i) { // get ith vertex
    return verts_[i];
}

Vec3f Model::vert(int f, int v) { // get ith vertex
    return verts_[faces_[f][v*3]];
}

Vec3f Model::tvert(int i) { // get ith texture vertex
    return tverts_[i];
}

Vec3f Model::tvert(int f,int v) { // get ith texture vertex
    return tverts_[faces_[f][v*3+1]];
}

Vec3f Model::nvert(int i) { // get ith normal
    return nverts_[i];
}
