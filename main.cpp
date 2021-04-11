#include "tgaimage.h"
#include <stdlib.h>

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

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
	float slope = std::abs(dy/(float)dx); // rise/run.
	float error = 0; // quantifies distance between y and the closest point on the straight line.
	for (int x=x0; x<x1; x++){
		if(std::abs(error+=slope)> .5){
			error-=1.;
			y += (y1 > y0 ? 1 : -1); //
		}
		steep ? image.set(y,x,color) : image.set(x,y,color);
	}
}

int main(int argc, char** argv) {
	TGAImage image(100, 100, TGAImage::RGB);

	//draw a line from (x0,y0) to (x1,y1)
	for (int i=0;i<1e7;i++){
		line(80,40,13,20, image, red);
		line(40,80,20,13, image, white);
	}
		image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
		image.write_tga_file("output.tga");
	return 0;
}
