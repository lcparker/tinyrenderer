#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"
#include "my_gl.h"

/* TODO list:
 * 	- performance improvements
 * 	- SIMD, OpenMP pragma
 * 	- remove duplicate methods (fill_shadow_buffer, triangle_blank)
 */

Model *model = NULL;

// WIDTH AND HEIGHT OF THE IMAGE BEING DRAWN
const int width  = 1000;
const int height = 1000;

// Depth buffer for camera and light reference frame - for hidden face removal
// and hard shadows, respectively.
float zbuffer[width*height] 	  {std::numeric_limits<float>::min()}; 
float shadow_buffer[width*height] {std::numeric_limits<float>::min()};

// POSITIONS: ORIGIN, CAMERA, CAMERA ORIENTATION, INCOMING LIGHT DIRECTION
const Vec3f centre(0,0,0);
const Vec3f eye(1,0,3);
const Vec3f up(0,1,0);
Vec3f light_dir = Vec3f(1,1,1);

// These messy variables speed up matrix multiplication in the core loop by quite a bit
// TODO see if it's quicker to pass by reference a couple times (it may be, depending on if these global vars are stored in memory
Matrix Minvtr; 
Matrix N;

struct Phong : public Shader {
	Vec3f verts[3];
	Vec3f dirs[3];
	Vec3f vert_textures[3];

	Vec3f vertex(int f, int v){
		/* Loads relevant shader data for the v'th vertex of face f. */
		Vec3f vert = model->vert(f,v);
		verts[v] = vert;
		dirs[v]  = eye-vert;
		vert_textures[v] = model->tvert(f,v);
		return Vec3f(N * Vec4f(vert)); // Make this more elegant
	}

	void fragment(Vec3f bary, TGAColor &c, Matrix &M, TGAImage &image){
		/* Interpolate the vertices according to the given barycentric
		 * coodinates, setting the correct colour to be drawn at that pixel. */
		Vec3f v    = verts[0]*bary[0]+verts[1]*bary[1]+verts[2]*bary[2];
		Vec3f pt_c = N * Vec4f(v);
		Vec3i pt_l = Vec3f(M*Vec4f(v)); // transform to light ray frame

		Vec3f n   = Minvtr*get_normal(vert_textures, bary, model->normals);
		Vec3f dir = (dirs[0]*bary[0] + dirs[1]*bary[1] + dirs[2]*bary[2]).normalize();
		Vec3f l   = ModelView*light_dir; 
	  	n.normalize();
		l.normalize();
		Vec3f r = (n*2.*(n*l)-l); // |R| = 1

		float diffuse  = std::max<float>(0.f, n*light_dir);
		float specular = std::max<float>(0.f,std::pow(dir*r,10));
		float shadow   = .3+.7*(shadow_buffer[pt_l.x*width+pt_l.y]<pt_l.z+50);

		float occlusion = 0.;	
		for(float a=0.; a<2*M_PI-1e-4; a+=M_PI/4){
			occlusion += M_PI/2-max_elevation_angle(zbuffer,
				   Vec2f(pt_c.x,pt_c.y),
				   Vec2f(cos(a),sin(a)),image);
		}
		occlusion /= (M_PI/2)*8;

		c = get_texture(vert_textures, bary, model->textures)*shadow*occlusion*(diffuse+0.6*specular)+5;
	}
};

int main(int argc, char** argv) {
	
	TGAImage image(width, height, TGAImage::RGB);

	std::string const& modelname = "obj/african_head";
	model = new Model(modelname);

	ViewPort = viewport(width/8., height/8., width*3/4., height*3/4., 255); 
	ModelView = view_frame(light_dir, centre, up);
	Projection = perspective((light_dir-centre).norm());
	Matrix world_to_light = ViewPort * Projection * ModelView; // Store coordinate transform

	// Fill the shadow buffer to determine shadow placement
	Phong bs;
	Vec3f pts[3];
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		for (int j=0; j<3; j++) pts[j] = bs.vertex(i,j);
		fill_shadow_buffer(pts, shadow_buffer, image);
	}

	ModelView = view_frame(eye, centre, up);
	Projection = perspective((eye-centre).norm());
	
	Minvtr = ModelView.inverse_transpose();
	N = ViewPort * Projection * ModelView;

	// Preset the z-buffer for ambient occlusion
	std::vector<int> face;
	for(int i=0; i<model->nfaces();i++){
 		face = model->face(i);
		for (int j=0; j<3; j++){
			pts[j] = bs.vertex(i,j);
		}	
		triangle_blank(pts, zbuffer, image);
	}

	light_dir.normalize();
	Phong shader;
	for(int i=0; i<model->nfaces();i++){
		Vec3f pts[3];
		for (int j=0; j<3; j++) pts[j] = shader.vertex(i,j);
		triangle(pts, shader, world_to_light, zbuffer, shadow_buffer, image);
	}

	delete model;

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}

