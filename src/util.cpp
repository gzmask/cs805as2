#include "util.h"
#include <math.h>

//pixel iterator for img panel.
ImagePanel foreach_pixel_exec(ImagePanel img, std::function<int(Ray)> ray_func){
  int i = 0;
  for (auto& pixel: img) { //foreach pixel in empty_img
    //using to_2d function to get x,y camera coordinates
    auto cam_xy = to_2d(i);

    //construct Ray
    Ray ray = ray_construction(cam_xy[0], cam_xy[1]);
    pixel = ray_func(ray);
    i++;
  }
  return img;
}

Ray ray_construction(int x, int y) {
  //transform VRP to world coordinate
    return {-1,-1,-1,
            -1,-1,-1};
}

Point mul(Matrix m, Point x) {
  return mul(x, m);
}

Point mul(Point x, Matrix m) {
  return {x[0]*m[0][0]+x[1]*m[0][1]+x[2]*m[0][2]+m[0][3],
          x[0]*m[1][0]+x[1]*m[1][1]+x[2]*m[1][2]+m[1][3],
          x[0]*m[2][0]+x[1]*m[2][1]+x[2]*m[2][2]+m[2][3]};
}



//initialize img panel to all 0s
ImagePanel init_img_panel(ImagePanel img) {
  for (auto& pixel: img) { //foreach pixel in empty_img
    pixel = 0;
  }
  return img;
}

//translate ray equation to an 0~255 shading value
int ray_tracing(Ray ray) {
  Intersection p = ray_objects_intersection(ray);
  return shading(p); 
}

//calculate the ray object intersection point
Intersection ray_objects_intersection(Ray ray) {
  auto sphere_hit = ray_sphere_intersection(ray);
  auto polygon_hit = ray_polygon_intersection(ray);
  if (sphere_hit.kd < 0 && polygon_hit.kd < 0) {
    return {-1,-1,-1,
            -1,-1,-1,
            -1.0};
  } else if (polygon_hit.kd < 0) {
    return sphere_hit; 
  } else if (sphere_hit.kd < 0) {
    return polygon_hit; 
  } else if (closer(sphere_hit.intersection, polygon_hit.intersection, {0,0,0})) {
    return sphere_hit; 
  } else {
    return polygon_hit; 
  }
}

Intersection ray_sphere_intersection(Ray ray) {
    return {-1,-1,-1,
            -1,-1,-1,
            -1.0};
}

Intersection ray_polygon_intersection(Ray ray) {
    return {-1,-1,-1,
            -1,-1,-1,
            -1.0};
}

//calculate shading value from 0~255 accordingly to intersection info
int shading(Intersection p) {
  if (p.kd < 0) {
    return -1;
  }

  return 255;
}

//==========helpers==========

//return if p1 is closer to p0 than p2
bool closer(Point p1, Point p2, Point p0) {
  float d1 = (p1[0] - p0[0])+(p1[1] - p0[1])+(p1[2] - p0[2]);
  float d2 = (p2[0] - p0[0])+(p2[1] - p0[1])+(p2[2] - p0[2]);
  return d1 < d2;
}


//Translate 2D array index of row column to 1D index.
//Notice that x, or column index, starts with 0. 
//If return value is -1 then there is an out-of-bounce error.
int to_1d(int x, int y) {
  if (x >= IMG_X || x < 0)
    return -1;
  if (y >= IMG_Y || y < 0)
    return -1;
  return (IMG_Y*y + x);
}

//Translate 1d array index to 2d
std::array<int, 2> to_2d(int x) {
  if (x>=(IMG_X*IMG_Y) || x < 0) {
    return {-1,-1};
  }
  int y_ = x / IMG_X; 
  int x_ = x % IMG_X;
  return {x_, y_};
}

//prints the img panel
void print_img_panel(ImagePanel img) {
  std::cout<<std::endl;
  for (auto& pixel : img) {
    std::cout<<pixel<<", ";
  }
  std::cout<<std::endl<<"Array size: "<<img.size()<<std::endl;
}

