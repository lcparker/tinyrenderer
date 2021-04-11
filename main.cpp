#include "tgaimage.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) {
	int x, y;
	// approximate line with a series of dots
	for (float t=0; t<1.0; t+=0.05){
		x = x0 + (x1-x0)*t;
		y = y0 + (y1-y0)*t;
		image.set(x, y, color); // draws a color dot at (x,y) 

	}
}

int main(int argc, char** argv) {
	TGAImage image(100, 100, TGAImage::RGB);

	//draw a line from (x0,y0) to (x1,y1)
	line(10,10,90,90, image, white);
	image.flip_vertically(); // i want to have the origin at the left bottom corner of the image
	image.write_tga_file("output.tga");
	return 0;
}
