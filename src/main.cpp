#include <iostream>
#include "util.h"

/* create a spherical object */
SPHERE obj1 = {1.0, 1.0, 1.0,	/* center of the circle */
		 1.0,		/* radius of the circle */
		 0.75};		/* diffuse reflection coefficient */

/* create a polygon object */
POLY4 obj2 = {	0.0, 0.0, 0.0,	/* v0 */
		0.0, 0.0, 2.0,	/* v1 */
		2.0, 0.0, 2.0,	/* v2 */
		2.0, 0.0, 0.0,	/* v3 */
		0.0, 1.0, 0.0,	/* normal of the polygon */
		0.8};		/* diffuse reflection coefficient */

//unsigned char img[ROWS][COLS];

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */
float xmin = 0.0175;
float ymin = -0.0175;
float xmax = -0.0175;
float ymax = 0.0175;


/* definition of the camera parameters */
Point VRP = {1.0, 2.0, 3.5};
Vector VPN = {0.0, -1.0, -2.5};
Vector VUP = {0.0, 1.0, 0.0};

float focal = 0.05;	/* focal length simulating 50 mm lens */

/* definition of light source */
Point LRP = {-10.0, 10.0, 2.0};	/* light position */
float Ip = 200.0;	/* intensity of the point light source */

/* Transformation from the world to the camera coordinates */
Matrix Mwc = get_M(VRP, VPN, VUP);
/* Transformation from the camera to the world coordinates */
Matrix Mcw = get_Mi(VRP, VPN, VUP);
/* Transformation from the world to light coordinates */
Matrix Mwl = get_T(LRP);
/* Transformation from the light to the world coordinates */
Matrix Mlw = get_Ti(LRP);

int main () {
  ImagePanel img;
  img = init_img_panel(img);
  img = foreach_pixel_exec(img, ray_tracing);
  //print_img_panel(img);

  //tests
  pmatrix("Mwc:", Mwc);

  pmatrix("Mcw:", Mcw);

  std::cout<<"to_1d function, expected to be 512:"<<std::endl;
  std::cout<<to_1d(0, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be 0, 1:"<<std::endl;
  std::cout<<to_2d(512)[0]<<std::endl;
  std::cout<<to_2d(512)[1]<<std::endl;

  std::cout<<"to_1d function, expected to be 513:"<<std::endl;
  std::cout<<to_1d(1, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be 1,1:"<<std::endl;
  std::cout<<to_2d(513)[0]<<std::endl;
  std::cout<<to_2d(513)[1]<<std::endl;

  std::cout<<"to_1d function, expected to be 1023:"<<std::endl;
  std::cout<<to_1d(511, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be 511,1:"<<std::endl;
  std::cout<<to_2d(1023)[0]<<std::endl;
  std::cout<<to_2d(1023)[1]<<std::endl;

  std::cout<<"to_1d function, expected to be -1:"<<std::endl;
  std::cout<<to_1d(512, 1)<<std::endl;

  std::cout<<"to_2d function, expected to be -1,-1:"<<std::endl;
  std::cout<<to_2d(512*512)[0]<<std::endl;
  std::cout<<to_2d(512*512)[1]<<std::endl;

  std::cout<<"closer function, expected to be 1 and 0:"<<std::endl;
  std::cout<<closer({1,1,1},{2,2,2},{0,0,0})<<std::endl;
  std::cout<<closer({3,3,3},{2,2,2},{0,0,0})<<std::endl;

  std::cout<<"mul function: point * matrix"<<std::endl;
  Point a = {3,3,3};
  std::cout<<mul(a, Mcw)[0]<<std::endl;
  std::cout<<mul(Mcw, a)[1]<<std::endl;
  std::cout<<mul(a, Mcw)[2]<<std::endl;

  pmatrix("mul funciton: matrix*matrix:", mul(Mwc, Mcw));

  std::cout<<"mul function: row * matrix"<<std::endl;
  Row b = {3,3,3,3};
  std::cout<<mul(b, Mcw)[0]<<std::endl;
  std::cout<<mul(Mcw, b)[1]<<std::endl;
  std::cout<<mul(b, Mcw)[2]<<std::endl;
  std::cout<<mul(b, Mcw)[3]<<std::endl;


  return 0;
}
