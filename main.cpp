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
const int width  = 1000;
const int height = 1000;

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

void triangle(Vec3f *v, float zbuffer[width][height], TGAImage &image, TGAColor color) {
	/* Using barycentric coordinates to draw the triangle, which
	 * allows us to exploit parallelisation, even though on a 
	 * single thread it will take longer. 
	 */

	// Upper left and lower right point of triangle	bounding box

	// whether these are int or float doesn't seem to affect performance.
	Vec2f ur = Vec2f(std::max({v[0].x,v[1].x,v[2].x}),std::max({v[0].y,v[1].y,v[2].y})); //upper right corner of bouding box
	Vec2f ll = Vec2f(std::min({v[0].x,v[1].x,v[2].x}),std::min({v[0].y,v[1].y,v[2].y})); //lower left corner

	/*
	for(int i=0;i<3;i++){
		std::cout << v[i] << std::endl;
	}
	*/
	for (int x = (ll.x+1.) * width/2.; x<=(ur.x+1.)*width/2.;x++){
		for (int y = (ll.y+1.)*height/2; y<=(ur.y+1.)*height/2.;y++){
			Vec2i pts[3] = {Vec2i((v[0].x+1.)*width/2,(v[0].y+1.)*height/2),
			   Vec2i((v[1].x+1.)*width/2,(v[1].y+1.)*height/2),
			   Vec2i((v[2].x+1.)*width/2,(v[2].y+1.)*height/2)};
			Vec3f comps  = barycentric(pts,Vec2i(x,y));
			if(comps.x < 0 || comps.y < 0 || comps.z < 0) continue;
			float z = (v[0].z)*comps.x + (v[1].z)*comps.y + (v[2].z)*comps.z;
			if(zbuffer[x][y] < z){
				zbuffer[x][y] = z;
		   		image.set(x,y,color);
			}
		}
	}
}

float brightness(Vec3f *v, Vec3f light_dirn){
	Vec3f normal = (v[2]-v[0])^(v[1]-v[0]);
	normal.normalize();
	return normal*light_dirn;
}

TGAColor get_texture(Vec2i p, Vec3f *tv, Vec3f bary, TGAImage &texture){
	int tw = texture.get_width();
	int th = texture.get_height();
	Vec2i txts[3] = {Vec2i(tv[0].x*tw,tv[0].y*th),
	   Vec2i(tv[1].x*tw,tv[1].y*th),
	   Vec2i(tv[2].x*tw,tv[2].y*th)};
	int tx = bary.x*txts[1].x + bary.y*txts[2].x + bary.z*txts[0].x;
	int ty = bary.x*txts[2].y + bary.y*txts[1].y + bary.z*txts[0].y;
	return texture.get(tx,ty);
}

void triangle_texture(Vec3f *v, Vec3f *tv, Vec3f light_dir, float zbuffer[width][height], TGAImage &image, TGAImage &texture) {

	Vec2f ur = Vec2f(std::max({v[0].x,v[1].x,v[2].x}),std::max({v[0].y,v[1].y,v[2].y})); //upper right corner of bouding box
	Vec2f ll = Vec2f(std::min({v[0].x,v[1].x,v[2].x}),std::min({v[0].y,v[1].y,v[2].y})); //lower left corner

	for (int x = (ll.x+1.) * width/2.; x<=(ur.x+1.)*width/2.;x++){
		for (int y = (ll.y+1.)*height/2; y<=(ur.y+1.)*height/2.;y++){
			Vec2i pts[3] = {Vec2i((v[0].x+1.)*width/2,(v[0].y+1.)*height/2),
			   				Vec2i((v[1].x+1.)*width/2,(v[1].y+1.)*height/2),
			   				Vec2i((v[2].x+1.)*width/2,(v[2].y+1.)*height/2)};
			Vec3f comps  = barycentric(pts,Vec2i(x,y));
			if(comps.x < 0 || comps.y < 0 || comps.z < 0) continue;
			float z = (v[0].z)*comps.x + (v[1].z)*comps.y + (v[2].z)*comps.z;
			if(zbuffer[x][y] < z){
				zbuffer[x][y] = z;
				float b = brightness(v, light_dir);
				if(b > 0){
					TGAColor color = get_texture(Vec2i(x,y), tv, comps, texture)*b;
					image.set(x,y,color);
				}
			}
		}
	}
}

void drawfacemesh(const char *modelname, const char *texturename, TGAImage &image){

	model = new Model(modelname); // face
	TGAImage texture;
    texture.read_tga_file(texturename); // The image containing the texture map
	texture.flip_vertically();

	float zbuffer[width][height];
	for(int i=0;i<width;i++){
		for(int j=0;j<height;j++){
			zbuffer[i][j] = -2;
		}
	}

	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		Vec3f pts[3];
		Vec3f txts[3]; // the three texture coordinates.
		for (int j=0; j<3; j++){
			pts[j] = (model->vert(face[j*3]));
			txts[j] = (model->tvert(face[(j*3)+1]));
		}	
		Vec3f light_dir = Vec3f(0,0,-1); // light goes in the -z dir'n
		triangle_texture(pts, txts, light_dir, zbuffer, image, texture);
	}

	delete model;
}

int main(int argc, char** argv) {
	
	TGAImage image(width, height, TGAImage::RGB);
	drawfacemesh("obj/african_head.obj", "obj/african_head_diffuse.tga", image);

	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}

