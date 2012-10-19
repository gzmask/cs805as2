#include <iostream>
#include <typeinfo>//debugging only
#include "util.h"

int main () {
  ImagePanel resultImg;
  resultImg = init_img_panel(resultImg);
  resultImg = foreach_pixel_exec(resultImg, [](int x){return x+2;});
  print_img_panel(resultImg);
  return 0;
}
