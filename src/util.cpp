#include "util.h"
#include <math.h>

ImagePanel foreach_pixel_exec(ImagePanel empty_img, std::function<int(int)> ray_func){
  for (auto& pixel: empty_img) { //foreach pixel in empty_img
    //std::cout<<"before: "<<pixel<<std::endl;
    pixel = ray_func(pixel);
    //std::cout<<"after: "<<pixel<<std::endl;
  }
  return empty_img;
}

ImagePanel init_img_panel(ImagePanel img) {
  for (auto& pixel: img) { //foreach pixel in empty_img
    pixel = 0;
  }
  return img;
}

//helpers
void print_img_panel(ImagePanel img) {
  std::cout<<std::endl;
  for (auto& pixel : img) {
    std::cout<<pixel<<", ";
  }
  std::cout<<std::endl<<"Array size: "<<img.size()<<std::endl;
}

