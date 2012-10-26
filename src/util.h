#ifndef UTIL_H
#define UTIL_H

//define preprocessing vars
#define IMG_X 512
#define IMG_Y 512
#define IMG_LEN ( IMG_X * IMG_Y )
/* definition of the image buffer */
#define ROWS IMG_Y
#define COLS IMG_X

#include <array>
#include <functional>
#include <iostream>
#include <iomanip>

//types
typedef std::array<int, IMG_LEN> ImagePanel;
typedef std::array<double, 3> Point;
typedef std::array<double, 3> Vector;
typedef std::array<Vector, 3> UVN;
typedef struct {
	Point intersection;	/* intersection point */
	Vector normal;	/* intersection polygon normal vector */
	double kd;	/* diffuse reflection coefficient of the surface */
} Intersection;
typedef struct {
	Point ref;	/* reference point, where the ray is from */
	Vector direction;	/* ray direction */
} Ray;
typedef std::array<double, 4> Row;
typedef std::array<Row, 4> Matrix;
typedef struct {
	double x, y, z;	/* center of the circle */
	double radius;	/* radius of the circle */
	double kd;	/* diffuse reflection coefficient */
} SPHERE;
typedef struct {
	double v[4][3];	/* list of vertices */
	double N[3];	/* normal of the polygon */
	double kd;	/* diffuse reflection coefficient */
} POLY4;

//functions
ImagePanel foreach_pixel_exec(ImagePanel, std::function<int(Ray)>);
ImagePanel init_img_panel(ImagePanel);
int ray_tracing(Ray);
Intersection ray_objects_intersection(Ray);
int shading(Intersection);
Intersection ray_sphere_intersection(Ray);
Intersection ray_polygon_intersection(Ray);
Ray ray_construction(int, int);

//helper functions
Point mul(Point, Matrix);
Point mul(Matrix, Point);
Matrix mul(Matrix, Matrix);
Row mul(Row, Matrix);
Row mul(Matrix, Row);
int to_1d(int, int);
std::array<int, 2> to_2d(int);
void print_img_panel(ImagePanel);
void pmatrix(std::string, Matrix);
bool closer(Point, Point, Point);
UVN get_uvn(Vector V1, Vector V2);
Matrix get_T(Point);
Matrix get_Ti(Point);
Matrix get_R(Point, Vector, Vector);
Matrix get_Ri(Point, Vector, Vector);
Matrix get_M(Point, Vector, Vector);
Matrix get_Mi(Point, Vector, Vector);
double get_length(Vector);
Vector cross_product(Vector, Vector);
Vector normalize(Vector);

//global vars
extern Matrix Mwc;
extern Matrix Rwc;
extern Matrix Twc;
extern Matrix Mcw;
extern Matrix Rcw;
extern Matrix Tcw;
extern Matrix Mwl;
extern Matrix Mlw;
extern double xmin;
extern double ymin;
extern double xmax;
extern double ymax;
extern Point VRP;
extern Vector VPN;
extern Vector VUP;
extern double focal;
extern Point LRP;
extern double Ip;
#endif
