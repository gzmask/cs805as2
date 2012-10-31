#include "util.h"
#include <cmath>

//pixel iterator for img panel.
ImagePanel foreach_pixel_exec(ImagePanel img, std::function<int(Ray)> ray_func) {
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

//ray constructor
Ray ray_construction(int x, int y) {
  //calculate x unit
  double x_delta = (xmax-xmin) / IMG_X;
  double y_delta = (ymax-ymin) / IMG_Y;

  //calculate the point on img panel with world coordinate
  double x_ = xmax - x_delta * x;
  double y_ = ymax - y_delta * y;

  //get vector v0. it is trival that VRP is p0
  Point p0 = VRP;
  Point p1_ = {x_, y_, focal};
  Point p1 = mul(Mcw, p1_);
  Vector v0_ = {p1[0] - p0[0],
    p1[1] - p0[1], 
    p1[2] - p0[2]};
  Vector v0 = normalize(v0_);

  /*
  if ((x==0 ) || (x==511)) {
    std::cout<<"img: x:"<<x<<", y:"<<y;
    //std::cout<<"p0: x:"<<VRP[0]<<", y:"<<VRP[1]<<", z:"<<VRP[2]<<std::endl;
    //std::cout<<"p0: x:"<<p0[0]<<", y:"<<p0[1]<<", z:"<<p0[2]<<std::endl;
    std::cout<<"p1_: x:"<<p1_[0]<<", y:"<<p1_[1]<<", z:"<<p1_[2]<<"=====";
    std::cout<<"p1: x:"<<p1[0]<<", y:"<<p1[1]<<", z:"<<p1[2]<<"=====";
    std::cout<<"v0: x:"<<v0[0]<<", y:"<<v0[1]<<", z:"<<v0[2]<<"=====";
    std::cout<<std::endl;
  }
  */

  return { p0[0], p0[1], p0[2],
           v0[0], v0[1], v0[2]};
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
  std::cout<<"Intersection: "<<p.intersection[0]<<","<<p.intersection[1]<<","<<p.intersection[2]<<"kd: "<<p.kd<<std::endl;
  return shading(p); 
}

//calculate the ray object intersection point
Intersection ray_objects_intersection(Ray ray) {
  auto sphere_hit = ray_sphere_intersection(ray, obj1);
  auto polygon_hit = ray_polygon_intersection(ray, obj2);
  if (sphere_hit.kd < 0 && polygon_hit.kd < 0) {
    return {-1,-1,-1,
            -1,-1,-1,
            -1.0};
  } else if (polygon_hit.kd < 0) {
    return sphere_hit; 
  } else if (sphere_hit.kd < 0) {
    return polygon_hit; 
  } else if (closer(sphere_hit.intersection, polygon_hit.intersection, ray.ref)) {
    return sphere_hit; 
  } else {
    return polygon_hit; 
  }
}

Intersection ray_sphere_intersection(Ray ray, SPHERE obj) {
  //get A,B,C
  //A = Xd^2 + Yd^2 + Zd^2
  double A = pow(ray.direction[0], 2) +
             pow(ray.direction[1], 2) +
             pow(ray.direction[2], 2);
  //B = 2 * (Xd * (X0 - Xc) + Yd * (Y0 - Yc) + Zd * (Z0 - Zc))
  double B = 2 * (ray.direction[0] * (ray.ref[0] - obj.x) +
                  ray.direction[1] * (ray.ref[1] - obj.y) +
                  ray.direction[2] * (ray.ref[2] - obj.z) );
  //C = (X0 - Xc)^2 + (Y0 - Yc)^2 + (Z0 - Zc)^2 - Sr^2
  double C = pow(ray.ref[0]-obj.x, 2) +
             pow(ray.ref[1]-obj.y, 2) +
             pow(ray.ref[2]-obj.z, 2) - 
             pow(obj.radius, 2);

  //get discriminant
  double discriminant = pow(B,2) - 4*C;

  //return null if discriminant is less than 0
  Intersection null_ = {-1,-1,-1,
                        -1,-1,-1,
                        -1.0};
  if (discriminant < 0)
    return null_;

  //compute t0 = (- B - (B^2 - 4*C)^1/2) / 2
  double t0 = (-B - sqrt(discriminant)) / 2;
  double t1 = (-B + sqrt(discriminant)) / 2;

  //compute the intersection point Ri = [x0 + xd * ti ,  y0 + yd * ti,  z0 + zd * ti]
  Point Ri;
  if (discriminant > 0) {
    Ri = {ray.ref[0] + ray.direction[0] * t0, 
          ray.ref[1] + ray.direction[1] * t0, 
          ray.ref[2] + ray.direction[2] * t0}; 
  } else {
    Ri = {ray.ref[0] + ray.direction[0] * t1, 
          ray.ref[1] + ray.direction[1] * t1, 
          ray.ref[2] + ray.direction[2] * t1}; 
  }

  //compute the surface normal SN = [(xi - xc)/Sr,   (yi - yc)/Sr,   (zi - zc)/Sr]
  Vector SN = { (Ri[0]-obj.x)/obj.radius,
                (Ri[1]-obj.y)/obj.radius,
                (Ri[2]-obj.z)/obj.radius };

  Intersection result = { Ri, SN, obj.kd };
  return result;
}

Intersection ray_polygon_intersection(Ray ray, POLY4 obj) {
  Intersection null_ = {-1,-1,-1, -1,-1,-1, -1.0};

  //compute ray plane intersection
  //first compute Pn · Rd = Vd 
  double Vd = dot_product(obj.N, ray.direction);

  //Vd = 0: ray is parallel to the plane
  if (Vd == 0) 
    return null_;

  //Vd > 0: plane facing away from the ray
  if (Vd > 0)
    return null_;

  //second compute V0 = -(Pn· R0 + D)
  double V0 = - (dot_product(obj.N, ray.ref) + 1);

  double t = V0/Vd;

  //If t < 0 then the ray intersects plane at the negative side of the ray
  if (t<0) 
    return null_;

  //compute intersection point: Pi = [Xi Yi Zi] = [X0 + Xd * t Y0 + Yd * t Z0 + Zd * t]
  Point Pi = { ray.ref[0] + ray.direction[0]*t, 
               ray.ref[1] + ray.direction[1]*t,
               ray.ref[2] + ray.direction[2]*t };

  //check if intersection point is inside the polygon
  if (in_poly4(Pi, obj))
    return { Pi, obj.N, obj.kd };
  else
    return null_;
}

//check if point is inside a 4 side polygon
bool in_poly4(Point p, POLY4 obj) {

  return true;
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
  //auto mt = get_T(vrp);

  //we can see vpn_ and vup_ as vectors. such that we can apply them to get_uvn function from q2
  auto uvn = get_uvn(vpn, vup);
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

//world to camera
Matrix get_M(Point vrp, Vector vpn, Vector vup) {
  return mul(get_R(vrp, vpn, vup), get_T(vrp));
}

//camera to world
Matrix get_Mi(Point vrp, Vector vpn, Vector vup) {
  return mul(get_Ti(vrp), get_Ri(vrp, vpn, vup));
}


//matrix multiplication
Point mul(Matrix m, Point x) {
  return mul(x, m);
}

Point mul(Point x, Matrix m) {
  double w =  m[3][0] * x[0]
        + m[3][1] * x[1]
        + m[3][2] * x[2]
        + m[3][3];
  return {(x[0]*m[0][0]+x[1]*m[0][1]+x[2]*m[0][2]+m[0][3])/w,
          (x[0]*m[1][0]+x[1]*m[1][1]+x[2]*m[1][2]+m[1][3])/w,
          (x[0]*m[2][0]+x[1]*m[2][1]+x[2]*m[2][2]+m[2][3])/w};
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
  Vector v1 = { (p1[0] - p0[0]), (p1[1] - p0[1]), (p1[2] - p0[2])};
  double d1 = get_length(v1);
  Vector v2 = { (p2[0] - p0[0]), (p2[1] - p0[1]), (p2[2] - p0[2])};
  double d2 = get_length(v2);
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

//dot product
double dot_product(Vector x, Vector y) {
  return x[0]*y[0]+x[1]*y[1]+x[2]*y[2]; 

}

//calculates cross product of two Vectors
Vector cross_product(Vector x, Vector y) {
  return { x[1]*y[2] - x[2]*y[1],
           x[2]*y[0] - x[0]*y[2],
           x[0]*y[1] - x[1]*y[0]};
}

//calculates length of a Vector
double get_length(Vector x) {
  return sqrt(pow(x[0],2)+pow(x[1],2)+pow(x[2],2));
}

