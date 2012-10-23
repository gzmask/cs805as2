#include <iostream>
#include <typeinfo>//debugging only
#include "util.h"

int main () {
  ImagePanel img;
  img = init_img_panel(img);
  img = foreach_pixel_exec(img, ray_tracing);
  //print_img_panel(img);

  //unit tests
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

  return 0;
}
