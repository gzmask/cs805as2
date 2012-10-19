#ifndef UTIL_H
#define UTIL_H

//define global vars
#define IMG_X 320
#define IMG_Y 240
#define IMG_LEN ( IMG_X * IMG_Y )

#include <array>
#include <functional>
#include <iostream>
typedef std::array<int, IMG_LEN> ImagePanel;
typedef std::array<float, 3> Ray;//assuming there are 3 parameters for ray equation

ImagePanel foreach_pixel_exec(ImagePanel, std::function<int(int)>);
ImagePanel init_img_panel(ImagePanel);

//helpers
void print_img_panel(ImagePanel);
#endif
