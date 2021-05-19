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
 * 	- SIMD, other code optimisations. Profiling
 * 	- access image, normals, textures from one TGAImage instance
 */

Model *model = NULL;

// WIDTH AND HEIGHT OF THE IMAGE BEING DRAWN
extern const int width  = 1000;
extern const int height = 1000;

float zbuffer[width*height] {std::numeric_limits<float>::min()}; // Might need to change this to -max if there's an error.
float shadow_buffer[width*height] {std::numeric_limits<float>::min()};

// POSITIONS: ORIGIN, CAMERA, CAMERA ORIENTATION, INCOMING LIGHT DIRECTION
const Vec3f centre(0,0,0);
const Vec3f eye(1,0,3);
const Vec3f up(0,1,0);
Vec3f light_dir = Vec3f(1,1,1);


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
		return Vec3f(ViewPort * Projection * ModelView * Vec4f(vert)); // Make this more elegant
	}

	bool fragment(Vec3f bary, TGAColor &c, Matrix &M, TGAImage &image){
		/* Interpolate the vertices according to the given barycentric
		 * coodinates, setting the correct colour to be drawn at that pixel.
		 * Returns false if vertex should not be drawn. */
		Vec3f v    = verts[0]*bary[0]+verts[1]*bary[1]+verts[2]*bary[2];
		Vec3f pt_c = ViewPort * Projection * ModelView * Vec4f(v);
		Vec3i pt_l = Vec3f(M*Vec4f(v)); // transform to light ray frame

		Vec3f n   = ModelView.inverse_transpose()*get_normal(vert_textures, bary, model->normals);
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
		return true;
	}
};

int main(int argc, char** argv) {
	
	TGAImage image(width, height, TGAImage::RGB);
	TGAImage sbuffer(width,height, TGAImage::GRAYSCALE);

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
		fill_shadow_buffer(pts, shadow_buffer, sbuffer);
	}

	ModelView = view_frame(eye, centre, up);
	Projection = perspective((eye-centre).norm());
	
	// Preset the z-buffer for ambient occlusion
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		Vec3f pts[3];
		for (int j=0; j<3; j++){
			pts[j] = bs.vertex(i,j);
		}	
		triangle_blank(pts, zbuffer, sbuffer);
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
	sbuffer.flip_vertically();
	sbuffer.write_tga_file("buffer.tga");
	return 0;
}

