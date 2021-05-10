#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green   = TGAColor(0, 255,   0,  255);
const TGAColor blue   = TGAColor(0, 0,   255,  255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
Vec3f eye(1,1,3);
Vec3f centre(0,0,0);
Vec3f up(0,1,0);
Vec3f light_dir = Vec3f(0,3,3); // light goes in the -z dir'n

// TODO: Get normal vectors from texturemap.

void line(Vec2i v0, Vec2i v1, TGAImage &image, TGAColor color) {
/* Draw a line by coloring the pixels between two points. */
	bool steep = false;
	if(std::abs(v0.x-v1.x) < std::abs(v0.y-v1.y)){
		std::swap(v0.x,v0.y);
		std::swap(v1.x,v1.y);
		steep = true; }
   	if(v1.x < v0.x) std::swap(v0,v1);

	// For speed, take as much stuff outside the loop as possible
	// Minimise number of operations
	int y = v0.y;
	int dy = v1.y-v0.y;
	int dx = v1.x-v0.x;
	// We can do this without floats - we're working with pixels, which have int coordinates.
	// Floating point arithmetic is slower (?) so we should avoid these if we can. 
	int dy2 = std::abs(dy)*2; 
	int yerror = 0;
	for (int x=v0.x; x<v1.x; x++){
		// yerror =k*dy for integer k. If k*dy>0.5dx, i.e., k*(dy/dx)>0.5 or k*2dy > dx,
		// then inc y and 'mod' the error
		if((yerror+=dy2) > dx){
			yerror-=2*dx;
			y += (v1.y > v0.y ? 1 : -1);
		}
		steep ? image.set(y,x,color) : image.set(x,y,color);
	}
}

void rasterize(Vec2i v0, Vec2i v1, int ybuffer[width], TGAImage &image, TGAColor color){
	// Draw a line, but only if it's closer to the camera than what's already
	// been drawn.
	if(v1.x < v0.x) std::swap(v0,v1);

	float t;
	for(int x = v0.x; x<v1.x;x++){
		t = (x-v0.x)/(float)(v1.x-v0.x);
		int y = (1.-t)*v0.y + t*v1.y;
		if (y > ybuffer[x]){
			ybuffer[x] = y;
			for(int i=0;i<30;i++) image.set(x,i,color);
		}
	}
}

Vec3f barycentric(Vec2i *pts, Vec2i P){
	Vec3f xs(pts[2].x-pts[0].x, pts[1].x-pts[0].x,pts[0].x-P.x);
	Vec3f ys(pts[2].y-pts[0].y, pts[1].y-pts[0].y,pts[0].y-P.y);
	Vec3f u = xs ^ ys;
	
	if(std::abs(u.z) <1) return Vec3f(-1,1,1); // triangle is "degenerate", just return negative
	return Vec3f(u.x/u.z,u.y/u.z,1.-(u.x+u.y)/u.z);
}

TGAColor get_texture(Vec3f *tv, Vec3f bary, TGAImage &texture){
	int tw = texture.get_width();
	int th = texture.get_height();
	Vec2i txts[3] = {Vec2i(tv[0].x*tw,tv[0].y*th),
	   Vec2i(tv[1].x*tw,tv[1].y*th),
	   Vec2i(tv[2].x*tw,tv[2].y*th)};
	int tx = bary.x*txts[2].x + bary.y*txts[1].x + bary.z*txts[0].x;
	int ty = bary.x*txts[2].y + bary.y*txts[1].y + bary.z*txts[0].y;
	return texture.get(tx,ty);
}

Vec3f get_normal(Vec3f *tv, Vec3f bary, TGAImage &normalmap){
	int nw = normalmap.get_width();
	int nh = normalmap.get_height();
	Vec2i normals[3] = {Vec2i(tv[0].x*nw,tv[0].y*nh),
	   					Vec2i(tv[1].x*nw,tv[1].y*nh),
						Vec2i(tv[2].x*nw,tv[2].y*nh)};
	int nx = bary.x*normals[2].x + bary.y*normals[1].x + bary.z*normals[0].x;
	int ny = bary.x*normals[2].y + bary.y*normals[1].y + bary.z*normals[0].y;
	TGAColor cnormal = normalmap.get(nx,ny); // normal encoded as RGB
	Vec3f normal = Vec3f((cnormal.r-128.)/256.,(cnormal.g-128.)/256.,(cnormal.b-128.)/256.);
	return normal.normalize();
}

void triangle(Vec3f *v, Vec3f *tv, Vec3f *dirs, Matrix &ModelView,Matrix &Ninv, Matrix &M, int zbuffer[width][height], int sbuffer[width][height], TGAImage &image, TGAImage &texture, TGAImage &normalmap) {

	Matrix ModelViewTrInv = (ModelView.transpose()).inverse();
	Vec3f Ml = ModelView * Matrix(light_dir);

	Vec2f ur = Vec2f(std::max({v[0].x,v[1].x,v[2].x}),std::max({v[0].y,v[1].y,v[2].y})); //upper right corner of bouding box
	Vec2f ll = Vec2f(std::min({v[0].x,v[1].x,v[2].x}),std::min({v[0].y,v[1].y,v[2].y})); //lower left corner

	// don't try to draw things that will be outside of the camera view
	for (int x = ll.x; x<=ur.x && x < width;x++){
		for (int y = ll.y; y<=ur.y && y < height;y++){
			Vec2i pts[3] = {Vec2i(v[0].x,v[0].y),
			   				Vec2i(v[1].x,v[1].y),
			   				Vec2i(v[2].x,v[2].y)};
			Vec3f comps  = barycentric(pts,Vec2i(x,y));
			if(comps.x < 0 || comps.y < 0 || comps.z < 0) continue;
			Vec3f this_v = v[0]*comps.z + v[1]*comps.y + v[2]*comps.x;
			if(zbuffer[x][y] < (int) this_v.z){
				Vec3f langle = M*Ninv*this_v;
					zbuffer[x][y] = (int) this_v.z;
					Vec3f normal = Vec3f(ModelViewTrInv * Matrix(get_normal(tv, comps, normalmap))).normalize();
					Vec3f dir = (dirs[0]*comps.z + dirs[1]*comps.y + dirs[2]*comps.x).normalize();
					Vec3f r = (normal*2.*(normal*Ml)-light_dir).normalize();
					float diffuse = std::max<float>(0.f, normal*Ml);
					//TODO: Get specular from a file
					float specular = std::max<float>(0.f,std::pow(dir*r,10));
			        float shadow = .3+.7*(sbuffer[(int)langle.x][(int)langle.y]<(int)langle.z+10); 
					TGAColor c = get_texture(tv, comps, texture);
					TGAColor colour;
					colour.r = std::min<float>(255,5 + c.r*shadow*(diffuse + 0.6*specular));
					colour.g = std::min<float>(255,5 + c.g*shadow*(diffuse + 0.6*specular));
					colour.b = std::min<float>(255,5 + c.b*shadow*(diffuse + 0.6*specular));

					image.set(x,y,colour);
			}
		}
	}
}

Matrix view_frame(Vec3f eye, Vec3f centre, Vec3f up){
	/* Generates the matrix which transforms vectors from the object frame to the 
	 * camera's frame */

	Vec3f z = (eye-centre).normalize();
	Vec3f x = (up^z).normalize();
	Vec3f y = (z^x).normalize();
	Matrix Minv  = Matrix::identity(4);
	Matrix Trans = Matrix::identity(4);
	for(int i = 0; i < 3; i++){
		Minv[0][i]  = x[i];
		Minv[1][i]  = y[i];
		Minv[2][i]  = z[i];
		Minv[i][3] = -centre[i];
	}

	return Minv*Trans;
}

Matrix perspective(float c){
	/* generate the perspective-view distortion matrix for v */

	Matrix result = Matrix::identity(4);	
	float camscale = -1./c;

	result[3][2] = camscale;

	return result;
}

Matrix viewport(int x, int y, int w, int h, int depth){

	Matrix result = Matrix::identity(4);
	result[0][0] = w/2.;
	result[1][1] = h/2.;
	result[2][2] = depth/2.;
	result[2][3] = depth/2.;
	result[0][3] = x + w/2.;
	result[1][3] = y + h/2.;

	return result;
}

void shadow_buffer(Vec3f *v, int sb[width][height], TGAImage &zbuf){
	Vec2f ur = Vec2f(std::max({v[0].x,v[1].x,v[2].x}),std::max({v[0].y,v[1].y,v[2].y})); //upper right corner of bouding box
	Vec2f ll = Vec2f(std::min({v[0].x,v[1].x,v[2].x}),std::min({v[0].y,v[1].y,v[2].y})); //lower left corner

	// don't try to draw things that will be outside of the camera view
	for (int x = ll.x; x<=ur.x && x < width;x++){
		for (int y = ll.y; y<=ur.y && y < height;y++){
			Vec2i pts[3] = {Vec2i(v[0].x,v[0].y),
			   				Vec2i(v[1].x,v[1].y),
			   				Vec2i(v[2].x,v[2].y)};
			Vec3f comps  = barycentric(pts,Vec2i(x,y));
			if(comps.x < 0 || comps.y < 0 || comps.z < 0) continue;
			float z = (v[0].z)*comps.z + (v[1].z)*comps.y + (v[2].z)*comps.x;
			if(sb[x][y] < (int) z){
				sb[x][y] = (int) z;
				TGAColor c;
				c.r = std::min<float>(255,(int)z);
				c.g = std::min<float>(255,(int)z);
				c.b = std::min<float>(255,(int)z);
				zbuf.set(x,y,c);
			}
		}
	}
}


void drawfacemesh(const char *modelname, const char *texturename, const char *normalsname, TGAImage &image){

	model = new Model(modelname); // face
	TGAImage texture;
	TGAImage normalmap;
    texture.read_tga_file(texturename); // The image containing the texture map
    normalmap.read_tga_file(normalsname); // The image containing the normals encoded in RGB coords
	texture.flip_vertically();
	normalmap.flip_vertically();

	int shadowbuffer[width][height] = {0};
	TGAImage zbuf(width, height, TGAImage::RGB);

	Matrix ModelView = view_frame(light_dir, centre, up);
	Matrix ViewPort = viewport(width/8., height/8., width*3/4., height*3/4., 255); // so the head won't take up the whole screen
	Matrix Projection = perspective((light_dir-centre).norm());
	Matrix M = ViewPort * Projection * ModelView;
	
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		Vec3f pts[3];
		for (int j=0; j<3; j++){
			Vec3f v = (model->vert(face[j*3]));
			Vec4f v4(v);
			pts[j] = Vec3f(M * Matrix(v4));
		}	
		shadow_buffer(pts, shadowbuffer, zbuf);
	}
	zbuf.flip_vertically();
	zbuf.write_tga_file("buffer.tga");

	int zbuffer[width][height] = {0}; // z-buffer in greyscale

	ModelView = view_frame(eye, centre, up);
	ViewPort = viewport(width/8., height/8., width*3/4., height*3/4., 255); // so the head won't take up the whole screen
	Projection = perspective((eye-centre).norm());

	Matrix Ninv = (ViewPort*Projection*ModelView).inverse();
	light_dir.normalize();

	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		Vec3f pts[3];
		Vec3f dirs[3]; // direction from vertex to camera
		Vec3f textures[3]; // the three texture coordinates.
		for (int j=0; j<3; j++){
			Vec3f v = (model->vert(face[j*3]));
			Vec4f v4(v);
			pts[j] = Vec3f(ViewPort * Projection * ModelView * Matrix(v4));
			dirs[j] = eye-v;
			textures[j] = (model->tvert(face[(j*3)+1]));
		}	
		triangle(pts, textures, dirs, ModelView, Ninv, M, zbuffer, shadowbuffer, image, texture, normalmap);
	}

	delete model;
}

int main(int argc, char** argv) {
	
	TGAImage image(width, height, TGAImage::RGB);
	drawfacemesh("obj/african_head.obj",
		   		 "obj/african_head_diffuse.tga",
				 "obj/african_head_nm.tga", image);

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}
