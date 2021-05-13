#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "tgaimage.h"
#include "geometry.h"
#include "model.h"

// TODO: Get normal vectors from texturemap.

Model *model = NULL;

// WIDTH AND HEIGHT OF THE IMAGE BEING DRAWN
const int width  = 1000;
const int height = 1000;

// POSITIONS: ORIGIN, CAMERA, CAMERA ORIENTATION, INCOMING LIGHT DIRECTION
Vec3f centre(0,0,0);
Vec3f eye(1.2,-0.8,3);
Vec3f up(0,1,0);
Vec3f light_dir = Vec3f(0,3,3);

Vec3f barycentric(Vec2i *pts, Vec2i P){
/* Finds the barycentric coordinates of a point P in a triangle with points
 * pts[i]. 
 * INPUT:
 * 	pts[i]: vertices of the triangle, projected onto the screen.
 * 	     P: a point for which the barycentric coordinates will be calculated.
 */

	Vec3f xs(pts[2].x-pts[0].x, pts[1].x-pts[0].x,pts[0].x-P.x);
	Vec3f ys(pts[2].y-pts[0].y, pts[1].y-pts[0].y,pts[0].y-P.y);
	Vec3f u = xs ^ ys;
	
	// If triangle is degenerate, just return negative, so triangle rasteriser
	// will not draw the pixel.
	if(std::abs(u.z) < .1) return Vec3f(-1,1,1); 
	return Vec3f(u.x/u.z,u.y/u.z,1.-(u.x+u.y)/u.z);
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
	Vec2i txts[3] = {Vec2i(tv[0].x*tw,tv[0].y*th),
	   				 Vec2i(tv[1].x*tw,tv[1].y*th),
	   				 Vec2i(tv[2].x*tw,tv[2].y*th)};
	int tx = bary.x*txts[2].x + bary.y*txts[1].x + bary.z*txts[0].x;
	int ty = bary.x*txts[2].y + bary.y*txts[1].y + bary.z*txts[0].y;
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
	Vec2i normals[3] = {Vec2i(nv[0].x*nw,nv[0].y*nh),
	   					Vec2i(nv[1].x*nw,nv[1].y*nh),
						Vec2i(nv[2].x*nw,nv[2].y*nh)};
	int nx = bary.x*normals[2].x + bary.y*normals[1].x + bary.z*normals[0].x;
	int ny = bary.x*normals[2].y + bary.y*normals[1].y + bary.z*normals[0].y;
	TGAColor cnormal = normalmap.get(nx,ny); // normal encoded as RGB
	Vec3f normal = Vec3f((cnormal.r-128.)/256.,(cnormal.g-128.)/256.,
						 (cnormal.b-128.)/256.);
	return normal.normalize();
}

float max_elevation_angle(float z[width][height], Vec2f pt, Vec2f dir){
/* Find the maximum angle of elevation from pt to another visible point in
 * direction dir. (For screen-based ambient occlusion.)
 * INPUT:
 *    z: the depth buffer.
 *   pt: the point we are finding the elevation angle from
 * 	dir: the direction in which we are looking for the elevation angle.
 */
 
	float max_slope = std::numeric_limits<float>::min();
	Vec2f d(0.,0.);
	for(float i = 1.f;i <100.;i+=1.){
		d = dir*i+pt;
		float run = dir.norm()*i;
		if(run < 1.) continue;
		if(d.x >=width || d.x <0 || d.y >= height || d.y < 0) return std::atan(max_slope);
		float rise = z[int(d.x)][int(d.y)]-z[int(pt.x)][int(pt.y)];
		max_slope = std::max(max_slope,rise/run);
	}
	
	return std::atan(max_slope);
}

void triangle(Vec3f *v, Vec3f *tv, Vec3f *dirs, Matrix &ModelView,Matrix &Ninv, Matrix &M, float zbuffer[width][height], float sbuffer[width][height], TGAImage &image, TGAImage &texture, TGAImage &normalmap) {
/* Rasterise all appropriate points inside a triangle specified by points v[i].
 * 
 * TODO: Clean and update.
 */
	//Matrix ModelViewTrInv = (ModelView.transpose()).inverse();
	//Vec3f Ml = ModelView * Matrix(light_dir);

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
			//if(zbuffer[x][y] <= (int) this_v.z){
			//	Vec3f langle = M*Ninv*this_v;
					//zbuffer[x][y] = (int) this_v.z;
					//Vec3f normal = Vec3f(ModelViewTrInv * Matrix(get_normal(tv, comps, normalmap))).normalize();
					//Vec3f dir = (dirs[0]*comps.z + dirs[1]*comps.y + dirs[2]*comps.x).normalize();
					//Vec3f r = (normal*2.*(normal*Ml)-light_dir).normalize();
					//float diffuse = std::max<float>(0.f, normal*Ml);
					//TODO: Get specular from a file
					//float specular = std::max<float>(0.f,std::pow(dir*r,10));
			        //float shadow = .3+.7*(sbuffer[(int)langle.x][(int)langle.y]<(int)langle.z+10); 
			//		TGAColor c = get_texture(tv, comps, texture);
//			}
		}
	}
}

void triangle_blank(Vec3f *v,float zbuffer[width][height], TGAImage &image) {
/* update the depth buffer as if the triangle specified by vertices v[i] were
 * to be drawn on the screen.
 *
 * TODO: Remove, packing this functionality into a shader.
 */

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
			float z= v[0].z*comps.z + v[1].z*comps.y + v[2].z*comps.x;
			if(zbuffer[x][y] < z){
				zbuffer[x][y] = z;
				image.set(x,y,TGAColor(z,z,z,1));
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

void shadow_buffer(Vec3f *v, float sb[width][height], TGAImage &zbuf){
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
			if(sb[x][y] < z){
				sb[x][y] = z;
				TGAColor c;
				c.r = std::min<float>(255,z);
				c.g = std::min<float>(255,z);
				c.b = std::min<float>(255,z);
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

	float shadowbuffer[width][height] = {-std::numeric_limits<float>::max()};
	TGAImage zbuf(width, height, TGAImage::RGB);

	Matrix ModelView = view_frame(light_dir, centre, up);
	Matrix ViewPort = viewport(width/8., height/8., width*3/4., height*3/4., 255); 
	Matrix Projection = perspective((light_dir-centre).norm());
	Matrix M = ViewPort * Projection * ModelView;

	
	/*
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
	*/

	float zbuffer[width][height] = {-std::numeric_limits<float>::max()}; // z-buffer in greyscale


	ModelView = view_frame(eye, centre, up);
	ViewPort = viewport(width/8., height/8., width*3/4., height*3/4., 255); 
	Projection = perspective((eye-centre).norm());

	Matrix Ninv = (ViewPort*Projection*ModelView).inverse();
	light_dir.normalize();
	
	// Need to preset z-buffer for ambient occlusion!
	for(int i=0; i<model->nfaces();i++){
		std::vector<int> face = model->face(i);
		Vec3f pts[3];
		for (int j=0; j<3; j++){
			Vec3f v = (model->vert(face[j*3]));
			Vec4f v4(v);
			pts[j] = Vec3f(ViewPort * Projection * ModelView * Matrix(v4));
		}	
		triangle_blank(pts, zbuffer,zbuf);
	}

	zbuf.flip_vertically();
	zbuf.write_tga_file("buffer.tga");

	/*
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
	*/

	// Occlusion step:
	float occlusion, a, inc;
	for(int x=0;x<width;x++){
		for(int y=0;y<height;y++){
				if(zbuffer[x][y] < 1e-5) continue;
				occlusion = 0.;
				for(a=0.; a<2*M_PI-1e-4; a+=M_PI/4){
					inc = M_PI/2 - max_elevation_angle(zbuffer, Vec2f(x,y), Vec2f(std::cos(a),std::sin(a)));
					occlusion += inc;
				}

				occlusion /= (M_PI/2)*8;
				if(occlusion < 0.005) std::cout << occlusion << std::endl;

				TGAColor colour;
				colour.r = std::min<float>(255,255*occlusion);
				colour.g = std::min<float>(255,255*occlusion);
				colour.b = std::min<float>(255,255*occlusion);

				image.set(x,y,TGAColor((occlusion*255),(occlusion*255),(occlusion*255),255));
		}

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
