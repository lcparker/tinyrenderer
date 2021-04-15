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
Model *model = NULL;
const int width  = 1000;
const int height = 1000;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
/* Draw a line by coloring the pixels between two points. */
	bool steep = false;
	if(std::abs(x0-x1) < std::abs(y0-y1)){
		std::swap(x0,y0);
		std::swap(x1,y1);
		steep = true; }
   	if(x1 < x0){ std::swap(x0,x1);
		std::swap(y0,y1);
	}

	// For speed, take as much stuff outside the loop as possible
	// Minimise number of operations
	int y = y0;
	int dy = y1-y0;
	int dx = x1-x0;
	// We can do this without floats - we're working with pixels, which have int coordinates.
	// Floating point arithmetic is slower (?) so we should avoid these if we can. 
	int dy2 = std::abs(dy)*2; 
	int yerror = 0;
	for (int x=x0; x<x1; x++){
		// yerror =k*dy for integer k. If k*dy>0.5dx, i.e., k*(dy/dx)>0.5 or k*2dy > dx,
		// then inc y and 'mod' the error
		if((yerror+=dy2) > dx){
			yerror-=2*dx;
			y += (y1 > y0 ? 1 : -1);
		}
		steep ? image.set(y,x,color) : image.set(x,y,color);
	}
}

Vec3f barycentric(Vec2i *pts, Vec2i P){
	Vec3f xs(pts[2].x-pts[0].x, pts[1].x-pts[0].x,pts[0].x-P.x);
	Vec3f ys(pts[2].y-pts[0].y, pts[1].y-pts[0].y,pts[0].y-P.y);
	Vec3f u = xs ^ ys;
	
	if(std::abs(u.z) <1) return Vec3f(-1,1,1); // triangle is "degenerate", just return negative
	return Vec3f(u.x/u.z,u.y/u.z,1.-(u.x+u.y)/u.z);
}

void triangle_bary(Vec2i *v, TGAImage &image, TGAColor color) {
	/* Using barycentric coordinates to draw the triangle, which
	 * allows us to exploit parallelisation, even though on a 
	 * single thread it will take longer. 
	 */

	// Upper left and lower right point of triangle	bounding box

	// whether these are int or float doesn't seem to affect performance.
	Vec2i ur = Vec2i(std::max({v[0].x,v[1].x,v[2].x}),std::max({v[0].y,v[1].y,v[2].y}));
	Vec2i ll = Vec2i(std::min({v[0].x,v[1].x,v[2].x}),std::min({v[0].y,v[1].y,v[2].y}));

	for (int x = ll.x; x<=ur.x;x++){
		for (int y = ll.y; y<=ur.y;y++){
			Vec3f comps  = barycentric(v,Vec2i(x,y));
			// neat trick: if cancels as soon as it can verify truth,
			// so it'll be faster to use || than &&
			if(comps.x < 0 || comps.y < 0 || comps.z < 0) continue;
		   	image.set(x,y,color);
		}
	}
}

void triangle_mono(Vec2i *v, TGAImage &image, TGAColor color) {
	// Draw a solid triangle. Uses an algorithm that's more suited
	// to single-thread processors.

	// About a 25% speedup by calling image.set() directly instead
	// of using line.
	if (v[1].x < v[0].x) std::swap(v[0],v[1]);
	if (v[2].x < v[0].x) std::swap(v[0],v[2]); 
	if (v[2].x < v[1].x) std::swap(v[1],v[2]); 

	int ya;
	int yb;
	
	std::cout << "drawing" << std::endl;

	for(int x = v[0].x; x<v[1].x;x++){
		float t = (x-v[0].x)/(float)(v[2].x-v[0].x);
		float s  = (x-v[0].x)/(float)(v[1].x-v[0].x);
		ya = (1.-t)*v[0].y + t*(v[2].y);		
		yb = (1.-s)*v[0].y + s*(v[1].y);
		std::cout << "x = " << x << std::endl;
		if(ya < yb) std::swap(ya,yb);	
		for(int y=yb;y<ya;y++){
		   	image.set(x,y,color);
		}
	}
	for(int x = v[2].x;x>=v[1].x;x--){
		float t = (v[2].x-x)/(float)(v[2].x-v[0].x);
		float s = (v[2].x-x)/(float)(v[2].x-v[1].x);
		ya = (1.-s)*v[2].y + s*v[1].y;		
		yb = (1.-t)*v[2].y + t*v[0].y;
		std::cout << "x = " << x << std::endl;
		if(ya < yb) std::swap(ya,yb);	
		for(int y=yb;y<ya;y++){
		   	image.set(x,y,color);
		}
	}
}

void drawfacemesh(const char *filename, TGAImage &image){
	model = new Model(filename);
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		Vec2i pts[3];
		Vec3f obj_coords[3]; // the coords as given in the obj file
		for (int j=0; j<3; j++){
			Vec3f v = model->vert(face[j]); 
			// obj file stores as "v a b c" where -1<= a,b,c<=1,
			// so need to rescale
			pts[j] = Vec2i((v.x+1.)*width/2.,(v.y+1.)*height/2.);
			obj_coords[j] = v;
		}	
		Vec3f normal = (obj_coords[2]-obj_coords[0])^(obj_coords[1]-obj_coords[0]);
		normal.normalize();
		Vec3f light_dirn = Vec3f(0,0,-1); // light goes in the -z dir'n
		float brightness = normal*light_dirn;
		if (brightness > 0){
			triangle_bary(pts,image,TGAColor(brightness*255,brightness*255,brightness*255,255));
		}
	}
	delete model;
}

int main(int argc, char** argv) {
	const char *filename;
	// the filename for the obj can optionally be passed as arg

	TGAImage image(width, height, TGAImage::RGB);

	if(argc==2){
		filename = argv[1];
	} else{
		filename = "obj/african_head.obj";
	}

	drawfacemesh(filename, image); // Draw a facemesh

	/*
	Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
	Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
	Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
	triangle_bary(t0, image, red); 
	triangle_bary(t1, image, white); 
	triangle_bary(t2, image, green);
	
	// this gives reasonable answers
	Vec3f b = barycentric(t0, Vec2i(20,90));

	std::cout << b << std::endl;

	*/
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}

