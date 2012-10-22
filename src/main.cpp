#include <iostream>
#include <typeinfo>//debugging only
#include "util.h"

int main () {
  ImagePanel resultImg;
  resultImg = init_img_panel(resultImg);
  resultImg = foreach_pixel_exec(resultImg, [](float x){return ray_tracing({x,x,x,x});});
  print_img_panel(resultImg);

  //unit tests
  std::cout<<to_1d(0, 1)<<std::endl;
  std::cout<<to_2d(512)[0]<<std::endl;
  std::cout<<to_2d(512)[1]<<std::endl;
  std::cout<<to_1d(1, 1)<<std::endl;
  std::cout<<to_2d(513)[0]<<std::endl;
  std::cout<<to_2d(513)[1]<<std::endl;
  std::cout<<to_1d(511, 1)<<std::endl;
  std::cout<<to_2d(1023)[0]<<std::endl;
  std::cout<<to_2d(1023)[1]<<std::endl;
  std::cout<<to_1d(512, 1)<<std::endl;
  std::cout<<to_2d(512*512)[0]<<std::endl;
  std::cout<<to_2d(512*512)[1]<<std::endl;
  return 0;
}
