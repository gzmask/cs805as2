#ifndef UTIL_H
#define UTIL_H

//define global vars
#define IMG_X 512
#define IMG_Y 512
#define IMG_LEN ( IMG_X * IMG_Y )

#include <array>
#include <functional>
#include <iostream>
typedef std::array<int, IMG_LEN> ImagePanel;
typedef std::array<float, 3> Point;
typedef std::array<float, 4> Ray;//assuming there are 3 parameters for ray equation

ImagePanel foreach_pixel_exec(ImagePanel, std::function<int(float)>);
ImagePanel init_img_panel(ImagePanel);
int ray_tracing(Ray);
Point ray_objects_intersection(Ray);

//helpers
int to_1d(int, int);
std::array<int, 2> to_2d(int);
void print_img_panel(ImagePanel);
#endif
