#include "bbox.h"

#include "GL/glew.h"

#include <algorithm>
#include <iostream>

namespace CGL {

double faceinter(const Ray& r, Vector3D n, Vector3D p) {
  return dot(p-r.o, n)/dot(r.d, n);
}

bool BBox::intersect(const Ray& r, double& t0, double& t1) const {

  // TODO (Part 2.2):
  // Implement ray - bounding box intersection test
  // If the ray intersected the bouding box within the range given by
  // t0, t1, update t0 and t1 with the new intersection times.
  Vector3D n[3];
  n[0] = Vector3D(1, 0, 0);
  n[1] = Vector3D(0, 1, 0);
  n[2] = Vector3D(0, 0, 1);
  double t[3][2], tmp;
  int i;
  double tmin = -1000000, tmax = 10000000000;
  for (i = 0 ; i < 3; i ++) {
    t[i][0] = faceinter(r, n[i], min);
    t[i][1] = faceinter(r, n[i], max);
    if (t[i][0]>t[i][1]) {
      tmp = t[i][0];
      t[i][0] = t[i][1];
      t[i][1] = tmp;
    }
    if (t[i][0]>tmin) {
      tmin = t[i][0];
    }
    if (t[i][1]<tmax)
      tmax = t[i][1];
  }
  

  if (tmin<=tmax) {
    t1 = tmax;
    t0 = tmin;
    return true;
  }

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
