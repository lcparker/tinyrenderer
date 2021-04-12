#include <vector>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"


const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;

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

void triangle(int x0, int y0, int x1, int y1, int x2, int y2, TGAImage &image, TGAColor color) {
	// draw a triangle
	line(x0,y0,x1,y1,image,color);
	line(x1,y1,x2,y2,image,color);
	line(x0,y0,x2,y2,image,color);
}

int main(int argc, char** argv) {
	if(argc==2){
		model = new Model(argv[1]);
	} else{
		model = new Model("obj/african_head.obj");
	}

	TGAImage image(width, height, TGAImage::RGB);
	// (copied) code to draw all the lines in the model
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i); //vector representing a given triangle
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
	
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	delete model;
	return 0;
}
