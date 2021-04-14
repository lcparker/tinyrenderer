#include <iostream>
#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
const TGAColor green   = TGAColor(0, 255,   0,  255);
Model *model = NULL;
const int width  = 180;
const int height = 180;

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

void triangle(Vec2i v0, Vec2i v1, Vec2i v2, TGAImage &image, TGAColor color) {
	// Draw a solid triangle.

	// About a 25% speedup by calling image.set() directly instead
	// of using line.
	if (v1.x < v0.x) std::swap(v0,v1);
	if (v2.x < v0.x) std::swap(v0,v2); 
	if (v2.x < v1.x) std::swap(v1,v2); 

	int ya = v0.y;
	int yb = v0.y;

	for(int x = v0.x; x<v1.x;x++){
		float t = (x-v0.x)/(float)(v2.x-v0.x);
		float s  = (x-v0.x)/(float)(v1.x-v0.x);
		ya = (1.-t)*v0.y + t*(v2.y);		
		yb = (1.-s)*v0.y + s*(v1.y);
		if(ya < yb) std::swap(ya,yb);	
		for(int y=yb;y<ya;y++){
		   	image.set(x,y,color);
		}
	}
	for(int x = v2.x;x>=v1.x;x--){
		float t = (v2.x-x)/(float)(v2.x-v0.x);
		float s = (v2.x-x)/(float)(v2.x-v1.x);
		ya = (1.-s)*v2.y + s*v1.y;		
		yb = (1.-t)*v2.y + t*v0.y;
		if(ya < yb) std::swap(ya,yb);	
		for(int y=yb;y<ya;y++){
		   	image.set(x,y,color);
		}
	}
}

void drawfacemesh(const char *filename, TGAImage &image){
// draw a face mesh in image according to model. (code copied)
	model = new Model(filename);
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		for (int j=0; j<3; j++){
			Vec3f v0 = model->vert(face[j]); 
			Vec3f v1 = model->vert(face[(j+1)%3]); 
			int x0 = (v0.x+1.)*width/2.;
			int y0 = (v0.y+1.)*height/2.; 
			int x1 = (v1.x+1.)*width/2.; 
			int y1 = (v1.y+1.)*height/2.; 
			line(x0, y0, x1, y1, image, white); 
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

	// drawfacemesh(filename, image); // Draw a facemesh

	Vec2i t0[3] = {Vec2i(10, 70),   Vec2i(50, 160),  Vec2i(70, 80)}; 
	Vec2i t1[3] = {Vec2i(180, 50),  Vec2i(150, 1),   Vec2i(70, 180)}; 
	Vec2i t2[3] = {Vec2i(180, 150), Vec2i(120, 160), Vec2i(130, 180)}; 
	line(t0[0].x,t0[0].y,t0[1].x,t0[1].y,image,white);
	line(t0[0].x,t0[0].y,t0[2].x,t0[2].y,image,white);
	line(t0[2].x,t0[2].y,t0[1].x,t0[1].y,image,white);
	line(t1[0].x,t1[0].y,t1[1].x,t1[1].y,image,white);
	line(t1[0].x,t1[0].y,t1[2].x,t1[2].y,image,white);
	line(t1[2].x,t1[2].y,t1[1].x,t1[1].y,image,white);
	line(t2[0].x,t2[0].y,t2[1].x,t2[1].y,image,white);
	line(t2[0].x,t2[0].y,t2[2].x,t2[2].y,image,white);
	line(t0[2].x,t0[2].y,t0[1].x,t0[1].y,image,white);
	for (int i=0; i<50000; i++){
		triangle(t0[0], t0[1], t0[2], image, red); 
		triangle(t1[0], t1[1], t1[2], image, white); 
		triangle(t2[0], t2[1], t2[2], image, green);
	}
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}
