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
  pmatrix("mwc:", Mwc);
    return {-1,-1,-1,
            -1,-1,-1};
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

//prints a matrix
void pmatrix(std::string str, Matrix m) {
  std::cout<<str<<std::endl;
  for (auto row : m) {
    for (auto num : row) {
      std::cout<<std::setw (10);
      std::cout<<num;
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

//get transformation matrix
Matrix get_T(Point vrp) {
  Row r1 = {1, 0, 0, -vrp[0]};
  Row r2 = {0, 1, 0, -vrp[1]};
  Row r3 = {0, 0, 1, -vrp[2]};
  Row r4 = {0, 0, 0, 1};
  return {r1, r2, r3, r4};
}

//get inverse transformation matrix
Matrix get_Ti(Point vrp) {
  Row r1 = {1, 0, 0, vrp[0]};
  Row r2 = {0, 1, 0, vrp[1]};
  Row r3 = {0, 0, 1, vrp[2]};
  Row r4 = {0, 0, 0, 1};
  return {r1, r2, r3, r4};
}

//get rotation matrix
Matrix get_R(Point vrp, Vector vpn, Vector vup) {
  //first get the translation matrix from world to view
  auto mt = get_T(vrp);

  //we can see vpn_ and vup_ as vectors. such that we can apply them to get_uvn function from q2
  auto uvn = get_uvn(vup, vpn);
  //finally contruct our roation matrix using method 2 on class notes
  Row r1 = { uvn[0][0],uvn[0][1],uvn[0][2],0 };
  Row r2 = { uvn[1][0],uvn[1][1],uvn[1][2],0 };
  Row r3 = { uvn[2][0],uvn[2][1],uvn[2][2],0 };
  Row r4 = { 0, 0, 0, 1 };
  return { r1, r2, r3, r4 };
}

//get inverse rotation matrix
Matrix get_Ri(Point vrp, Vector vpn, Vector vup) {
  Matrix m = get_R(vrp, vpn, vup);
  Row r1 = { m[0][0], m[1][0], m[2][0], m[3][0] };
  Row r2 = { m[0][1], m[1][1], m[2][1], m[3][1] };
  Row r3 = { m[0][2], m[1][2], m[2][2], m[3][2] };
  Row r4 = { m[0][3], m[1][3], m[2][3], m[3][3] };
  return {r1,r2,r3,r4};
}

Matrix get_M(Point vrp, Vector vpn, Vector vup) {
  return mul(get_R(vrp, vpn, vup), get_T(vrp));
}

Matrix get_Mi(Point vrp, Vector vpn, Vector vup) {
  return mul(get_Ti(vrp), get_Ri(vrp, vpn, vup));
}


//matrix multiplication
Point mul(Matrix m, Point x) {
  return mul(x, m);
}

Point mul(Point x, Matrix m) {
  return {x[0]*m[0][0]+x[1]*m[0][1]+x[2]*m[0][2]+m[0][3],
          x[0]*m[1][0]+x[1]*m[1][1]+x[2]*m[1][2]+m[1][3],
          x[0]*m[2][0]+x[1]*m[2][1]+x[2]*m[2][2]+m[2][3]};
}

Matrix mul(Matrix m, Matrix n) {
  Row r1 = {m[0][0]*n[0][0]+m[0][1]*n[1][0]+m[0][2]*n[2][0]+m[0][3]*n[3][0],
            m[0][0]*n[0][1]+m[0][1]*n[1][1]+m[0][2]*n[2][1]+m[0][3]*n[3][1],
            m[0][0]*n[0][2]+m[0][1]*n[1][2]+m[0][2]*n[2][2]+m[0][3]*n[3][2],
            m[0][0]*n[0][3]+m[0][1]*n[1][3]+m[0][2]*n[2][3]+m[0][3]*n[3][3]};
  Row r2 = {m[1][0]*n[0][0]+m[1][1]*n[1][0]+m[1][2]*n[2][0]+m[1][3]*n[3][0],
            m[1][0]*n[0][1]+m[1][1]*n[1][1]+m[1][2]*n[2][1]+m[1][3]*n[3][1],
            m[1][0]*n[0][2]+m[1][1]*n[1][2]+m[1][2]*n[2][2]+m[1][3]*n[3][2],
            m[1][0]*n[0][3]+m[1][1]*n[1][3]+m[1][2]*n[2][3]+m[1][3]*n[3][3]};
  Row r3 = {m[2][0]*n[0][0]+m[2][1]*n[1][0]+m[2][2]*n[2][0]+m[2][3]*n[3][0],
            m[2][0]*n[0][1]+m[2][1]*n[1][1]+m[2][2]*n[2][1]+m[2][3]*n[3][1],
            m[2][0]*n[0][2]+m[2][1]*n[1][2]+m[2][2]*n[2][2]+m[2][3]*n[3][2],
            m[2][0]*n[0][3]+m[2][1]*n[1][3]+m[2][2]*n[2][3]+m[2][3]*n[3][3]};
  Row r4 = {m[3][0]*n[0][0]+m[3][1]*n[1][0]+m[3][2]*n[2][0]+m[3][3]*n[3][0],
            m[3][0]*n[0][1]+m[3][1]*n[1][1]+m[3][2]*n[2][1]+m[3][3]*n[3][1],
            m[3][0]*n[0][2]+m[3][1]*n[1][2]+m[3][2]*n[2][2]+m[3][3]*n[3][2],
            m[3][0]*n[0][3]+m[3][1]*n[1][3]+m[3][2]*n[2][3]+m[3][3]*n[3][3]};
  return {r1,r2,r3,r4};
}

Row mul(Row x, Matrix m) {
  return {x[0]*m[0][0]+x[1]*m[0][1]+x[2]*m[0][2]+x[3]*m[0][3],
          x[0]*m[1][0]+x[1]*m[1][1]+x[2]*m[1][2]+x[3]*m[1][3],
          x[0]*m[2][0]+x[1]*m[2][1]+x[2]*m[2][2]+x[3]*m[2][3],
          x[0]*m[3][0]+x[1]*m[3][1]+x[2]*m[3][2]+x[3]*m[3][3]};
}

Row mul(Matrix m, Row x) {
  return mul(x, m);
}



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

//get u,v,n from two non-collinear vectors
UVN get_uvn(Vector V1, Vector V2) {

  //get n, which is just normalized V1
  Vector n = normalize(V1); 

  //get u, which is normalized V2 x V1
  Vector u = normalize(cross_product(V2, V1));

  //get v, which is normalized n x u
  Vector v = normalize(cross_product(n, u));

  return {u,v,n};
}

//normalize a Vector
Vector normalize(Vector x) {
  return { x[0]/get_length(x), 
           x[1]/get_length(x), 
           x[2]/get_length(x) }; 
}

//calculates cross product of two Vectors
Vector cross_product(Vector x, Vector y) {
  return { x[1]*y[2] - x[2]*y[1],
           x[2]*y[0] - x[0]*y[2],
           x[0]*y[1] - x[1]*y[0]};
}

//calculates length of a Vector
float get_length(Vector x) {
  return sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
}

