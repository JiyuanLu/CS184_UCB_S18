#include <nanogui/nanogui.h>

#include "../clothMesh.h"
#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(PointMass &pm) {
  // TODO (Part 3.1): Handle collisions with spheres.
	Vector3D direction = pm.position - origin;
	double length = direction.norm();
	if(length <= radius){
		direction.normalize();
		Vector3D tangent_direction = direction * radius + origin - pm.last_position;
		pm.position = pm.last_position + tangent_direction * (1.0 - friction);
	}
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  Misc::draw_sphere(shader, origin, radius * 0.92);
}
