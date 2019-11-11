#include "sphere.h"

#include <cmath>

#include  "../bsdf.h"
#include "../misc/sphere_drawing.h"

namespace CGL { namespace StaticScene {

bool Sphere::test(const Ray& r, double& t1, double& t2) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  double a = dot(r.d, r.d);
  double b = 2*dot(r.o-o, r.d);
  double c = dot(r.o-o, r.o-o)-r2;
  double delta = b*b-4*a*c;
  if (delta>=0) {
    t1 = (-b-sqrt(delta))/(2*a);
    t2 = (-b+sqrt(delta))/(2*a);
    if (t1>=r.min_t&&t1<=r.max_t) {
      r.max_t = t1;
      return true;
    }
    
    if (t2<=r.max_t && t2>= r.min_t && t1 <=0) {
      r.max_t = t2;
      return true;
    }
    
  }
  return false;
}

bool Sphere::intersect(const Ray& r) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note that you might want to use the the Sphere::test helper here.
  double t1, t2;
  bool f = test(r, t1, t2);
  if (f) 
    return true;
  return false;

}

bool Sphere::intersect(const Ray& r, Intersection *i) const {

  // TODO (Part 1.4):
  // Implement ray - sphere intersection.
  // Note again that you might want to use the the Sphere::test helper here.
  // When an intersection takes place, the Intersection data should be updated
  // correspondingly.
  double t1, t2;
  bool f = test(r, t1, t2);
  if (f) {
    if (t1>0)
      i->t = t1;
    else i->t = t2;
    Vector3D p = r.o+i->t*r.d;
    Vector3D n = p-o;
    n.normalize();
    i->n = n;
    i->primitive = this;
    i->bsdf = get_bsdf();
    return true;
  }
  return false;

}

void Sphere::draw(const Color& c) const {
  Misc::draw_sphere_opengl(o, r, c);
}

void Sphere::drawOutline(const Color& c) const {
    //Misc::draw_sphere_opengl(o, r, c);
}


} // namespace StaticScene
} // namespace CGL
