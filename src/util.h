#ifndef UTIL_H
#define UTIL_H

//define preprocessing vars
#define IMG_X 512
#define IMG_Y 512
#define IMG_LEN ( IMG_X * IMG_Y )

#include <array>
#include <functional>
#include <iostream>
typedef std::array<int, IMG_LEN> ImagePanel;
typedef std::array<float, 3> Point;
typedef std::array<float, 3> Vector;
typedef struct {
	Point intersection;	/* intersection point */
	Vector normal;	/* intersection polygon normal vector */
	float kd;	/* diffuse reflection coefficient of the surface */
} Intersection;
typedef struct {
	Point ref;	/* reference point, where the ray is from */
	Vector direction;	/* ray direction */
} Ray;
typedef std::array<float, 4> Row;
typedef std::array<Row, 4> Matrix;

ImagePanel foreach_pixel_exec(ImagePanel, std::function<int(Ray)>);
ImagePanel init_img_panel(ImagePanel);
int ray_tracing(Ray);
Intersection ray_objects_intersection(Ray);
int shading(Intersection);
Intersection ray_sphere_intersection(Ray);
Intersection ray_polygon_intersection(Ray);
Ray ray_construction(int, int);

//helpers
Point mul(Point, Matrix);
Point mul(Matrix, Point);
int to_1d(int, int);
std::array<int, 2> to_2d(int);
void print_img_panel(ImagePanel);
bool closer(Point, Point, Point);
#endif
