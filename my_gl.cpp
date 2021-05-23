#include <iostream>
#include <vector>
#include <cmath>
//#include <omp.h>
#include "my_gl.h"

Matrix ModelView;
Matrix ViewPort;
Matrix Projection;

Vec3f barycentric(Vec2f *pts, Vec2f P){
/* Finds the barycentric coordinates of a point P in a triangle with points
 * pts[i]. 
 * INPUT:
 * 	pts[i]: vertices of the triangle, projected onto the screen.
 * 	     P: a point for which the barycentric coordinates will be calculated.
 */

	Vec3f xs(pts[2][0]-pts[0][0], pts[1][0]-pts[0][0],pts[0][0]-P[0]);
	Vec3f ys(pts[2][1]-pts[0][1], pts[1][1]-pts[0][1],pts[0][1]-P[1]);
	Vec3f u = xs ^ ys;
	
	// If triangle is degenerate, just return negative, so triangle rasteriser
	// will not draw the pixel.
	if(std::abs(u.z) < .1) return Vec3f(-1,1,1); 
	return Vec3f(1.-(u.x+u.y)/u.z,u.y/u.z,u.x/u.z);
}

TGAColor get_texture(Vec3f *tv, Vec3f bary, TGAImage &texture){
/* Get the texture for the pixel to be drawn.
 * INPUT:
 *    tv[i]: vertices, to be converted to coordinates in the texture map.
 *     bary: the barycentric coordinates of the point to be drawn.
 *  texture: the texturemap image containing the textures for the object.
 */
	int tw = texture.get_width();
	int th = texture.get_height();
	// Pixel coordinates of vertices tv[i]
	Vec2i txts[3] = {Vec2i(tv[0][0]*tw,tv[0][1]*th),
	   				 Vec2i(tv[1][0]*tw,tv[1][1]*th),
	   				 Vec2i(tv[2][0]*tw,tv[2][1]*th)};
	int tx = bary[0]*txts[0][0] + bary[1]*txts[1][0] + bary[2]*txts[2][0];
	int ty = bary[0]*txts[0][1] + bary[1]*txts[1][1] + bary[2]*txts[2][1];
	return texture.get(tx,ty);
}

Vec3f get_normal(Vec3f *nv, Vec3f bary, TGAImage &normalmap){
/* Finds the normal for the point specified by the given barycentric
 * coordinates.
 * INPUT:
 *      nv[i]: vertices to be converted to coordinates in the normal map.
 *       bary: the barycentric coordinates of the point to be drawn.
 *  normalmap: image containing the normals for the object.
 */
	int nw = normalmap.get_width();
	int nh = normalmap.get_height();
	// Pixel coordinate form of nv[i]
	Vec2i normals[3] = {Vec2i(nv[0][0]*nw,nv[0][1]*nh),
	   					Vec2i(nv[1][0]*nw,nv[1][1]*nh),
						Vec2i(nv[2][0]*nw,nv[2][1]*nh)};
	int nx = bary[0]*normals[0][0] + bary[1]*normals[1][0] + bary[2]*normals[2][0];
	int ny = bary[0]*normals[0][1] + bary[1]*normals[1][1] + bary[2]*normals[2][1];
	TGAColor cnormal = normalmap.get(nx,ny); // normal encoded as RGB
	Vec3f normal = Vec3f((cnormal.r-128.)/256.,(cnormal.g-128.)/256.,
						 (cnormal.b-128.)/256.);
	return normal.normalize();
}

//TODO: remove image dependency, only needed to get height
float max_elevation_angle(float *z, Vec2f pt, Vec2f dir, TGAImage &image){
/* Find the maximum angle of elevation from pt to another visible point in
 * direction dir. (For screen-based ambient occlusion.)
 * INPUT:
 *    z: the depth buffer.
 *   pt: the point we are finding the elevation angle from
 * 	dir: the direction in which we are looking for the elevation angle.
 */
 
	int w = image.get_width();
	int h = image.get_height();
	float max_slope = std::numeric_limits<float>::min();
	Vec2f d(0.,0.);
	for(float i = 1.f;i <5.;i+=1.){
		d = dir*i+pt;
		float run = dir.norm()*i;
		if(d.x >= w|| d.x <0 || d.y >= h || d.y < 0) return std::atan(max_slope);
		float rise = z[int(d.x)*w+int(d.y)]-z[int(pt.x)*w+int(pt.y)];
		max_slope = std::max(max_slope,rise/run);
	}
	
	return std::atan(max_slope);
}

void triangle(Vec3f v[3], Shader &shader, Matrix &world_to_light , float *zbuffer, float *sbuffer, TGAImage &image) {
/* Rasterise all appropriate points inside a triangle specified by points v[i].
 * 
 * TODO: Clean and update.
 */
	
	// Upper right and lower left corners of bounding box
	Vec2f ur = Vec2f(std::max({v[0][0],v[1][0],v[2][0]}),
					 std::max({v[0][1],v[1][1],v[2][1]})); 
	Vec2f ll = Vec2f(std::min({v[0][0],v[1][0],v[2][0]}),
					 std::min({v[0][1],v[1][1],v[2][1]})); 

	Vec2f pts[3] = {Vec2f(v[0][0],v[0][1]),Vec2f(v[1][0],v[1][1]),
					Vec2f(v[2][0],v[2][1])};
	// don't try to draw things that will be outside of the camera view
	int xmax = std::min<int>(ur.x,image.get_width());
	int ymax = std::min<int>(ur.y,image.get_height());
//#pragma omp parallel for collapse(2)
	for (int x = ll.x; x<=xmax;x++){
		for (int y = ll.y; y<=ymax;y++){
			Vec3f bary  = barycentric(pts,Vec2f(x,y));
			if(bary.x < 0 || bary.y < 0 || bary.z < 0) continue;
			Vec3f this_v = v[0]*bary[0] + v[1]*bary[1] + v[2]*bary[2];
			if(zbuffer[x*image.get_width()+y] <= this_v.z){
				TGAColor c;
				shader.fragment(bary,c,world_to_light,image);
				image.set(x,y,c);
			}
		}
	}
}

void triangle_blank(Vec3f *v,float *zbuffer, TGAImage &image) {
/* update the depth buffer as if the triangle specified by vertices v[i] were
 * to be drawn on the screen.
 *
 * TODO: Remove, packing this functionality into a shader.
 */

	// Upper right and lower left corners of the bounding box for triangle
	Vec2f ur = Vec2f(std::max({v[0][0],v[1][0],v[2][0]}),
					 std::max({v[0][1],v[1][1],v[2][1]}));
	Vec2f ll = Vec2f(std::min({v[0][0],v[1][0],v[2][0]}),
					 std::min({v[0][1],v[1][1],v[2][1]})); 

	Vec2f pts[3] = {Vec2f(v[0][0],v[0][1]),Vec2f(v[1][0],v[1][1]),
	   				Vec2f(v[2][0],v[2][1])};
	int w = image.get_width();
	for (int x = ll.x; x<=ur.x && x < image.get_width();x++){
		for (int y = ll.y; y<=ur.y && y < image.get_height();y++){
			Vec3f bary  = barycentric(pts,Vec2f(x,y));
			if(bary.x < 0 || bary.y < 0 || bary.z < 0) continue;
			float z= v[0][2]*bary[0] + v[1][2]*bary[1] + v[2][2]*bary[2];
			if(zbuffer[x*w+y] < z){
				zbuffer[x*w+y] = z;
				//image.set(x,y,TGAColor(z,z,z,1));
			}
		}
	}
}

Matrix view_frame(Vec3f this_eye, Vec3f this_centre, Vec3f this_up){
	/* Generates the matrix which transforms vectors from the object frame to the 
	 * camera's frame */

	Vec3f z = (this_eye-this_centre).normalize();
	Vec3f x = (this_up^z).normalize();
	Vec3f y = (z^x).normalize();
	Matrix Minv  = Matrix::identity(4);
	Matrix Trans = Matrix::identity(4);
	for(int i = 0; i < 3; i++){
		Minv[0][i]  = x[i];
		Minv[1][i]  = y[i];
		Minv[2][i]  = z[i];
		Minv[i][3] = -this_centre[i];
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

void fill_shadow_buffer(Vec3f *v, float *sb, TGAImage &buffer){
	Vec2f ur = Vec2f(std::max({v[0][0],v[1][0],v[2][0]}),std::max({v[0][1],v[1][1],v[2][1]})); //upper right corner of bouding box
	Vec2f ll = Vec2f(std::min({v[0][0],v[1][0],v[2][0]}),std::min({v[0][1],v[1][1],v[2][1]})); //lower left corner
	Vec2f pts[3] = {Vec2f(v[0][0],v[0][1]),Vec2f(v[1][0],v[1][1]),
					Vec2f(v[2][0],v[2][1])};

	// don't try to draw things that will be outside of the camera view
	int w = buffer.get_width();
	int h = buffer.get_height();
	for (int x = ll.x; x<=ur.x && x < w;x++){
		for (int y = ll.y; y<=ur.y && y < h;y++){
			Vec3f bary  = barycentric(pts,Vec2f(x,y));
			if(bary.x < 0 || bary.y < 0 || bary.z < 0) continue;
			float z = v[0][2]*bary[0] + v[1][2]*bary[1] + v[2][2]*bary[2];
			if(sb[x*w+y] < z){
				sb[x*w+y] = z;
				buffer.set(x,y,TGAColor(z,z,z,1));
			}
		}
	}
}

