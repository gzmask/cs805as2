#include "util.h"
#include <math.h>

ImagePanel foreach_pixel_exec(ImagePanel empty_img, std::function<int(float)> ray_func){
  for (auto& pixel: empty_img) { //foreach pixel in empty_img
    //std::cout<<"before: "<<pixel<<std::endl;
    pixel = ray_func(1.0);
    //std::cout<<"after: "<<pixel<<std::endl;
  }
  return empty_img;
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
  Point p = ray_objects_intersection(ray);
  return ray[0]+ray[1]+ray[2]+ray[3];
}

//calculate the ray object intersection point
Point ray_objects_intersection(Ray ray) {
  return {1,2,3};
}

//==========helpers==========

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

