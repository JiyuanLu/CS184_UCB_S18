#include "bbox.h"
// #include "math"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // Part 2, Task 2:
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.

  double tmin, tmax, tymin, tymax, tzmin, tzmax; 

  tmin = (min.x - r.o.x) / r.d.x; 
  tmax = (max.x - r.o.x) / r.d.x; 
  if (tmin > tmax)
    std::swap(tmin,tmax);
  tymin = (min.y - r.o.y) / r.d.y; 
  tymax = (max.y - r.o.y) / r.d.y; 
  if (tymin > tymax)
    std::swap(tymin,tymax);
  tzmin = (min.z - r.o.z) / r.d.z; 
  tzmax = (max.z - r.o.z) / r.d.z; 
  if (tzmin > tzmax)
    std::swap(tzmin,tzmax);

  tmin = std::max(std::max(tmin, tymin), tzmin);
  tmax = std::min(std::min(tmax, tymax), tzmax);
  t0 = tmin;
  t1 = tmax;
  if(tmin <= tmax && tmax >= 0)
    return true;
  else
    return false;
}

void BBox::draw(Color c) const {

  glColor4f(c.r, c.g, c.b, c.a);

  // top
  glBegin(GL_LINE_STRIP);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(max.x, max.y, max.z);
  glEnd();

  // bottom
  glBegin(GL_LINE_STRIP);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, min.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glEnd();

  // side
  glBegin(GL_LINES);
  glVertex3d(max.x, max.y, max.z);
  glVertex3d(max.x, min.y, max.z);
  glVertex3d(max.x, max.y, min.z);
  glVertex3d(max.x, min.y, min.z);
  glVertex3d(min.x, max.y, min.z);
  glVertex3d(min.x, min.y, min.z);
  glVertex3d(min.x, max.y, max.z);
  glVertex3d(min.x, min.y, max.z);
  glEnd();

}

std::ostream& operator<<(std::ostream& os, const BBox& b) {
  return os << "BBOX(" << b.min << ", " << b.max << ")";
}

} // namespace CGL
