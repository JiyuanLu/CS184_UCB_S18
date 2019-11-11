#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

//check if a point mass position is valid
bool checkValid(int x, int y, int num_width_points, int num_height_points){
  if(x >= 0 && x < num_width_points && y >=0 && y < num_height_points)
    return true;
  return false;
}

void Cloth::buildGrid() {
  // TODO (Part 1.1): Build a grid of masses.
  double unit_width = width / num_width_points;
  double unit_height = height / num_width_points;

  for(int j = 0; j < num_height_points; ++j){
    for(int i = 0; i < num_width_points; ++i){

      Vector3D position;
      bool pin;

      if(orientation == 0){  // horizontal
        position = Vector3D(unit_width * i, 1.0, unit_height * j);
      }
      else{                 // vertical
        double z = (rand() / (double) RAND_MAX) * 0.002 - 0.001;
        position = Vector3D(unit_width * i, unit_height * j, z);
      }

      if(std::find(pinned.begin(), pinned.end(), std::vector<int>{i, j}) != pinned.end()){
        pin = true;
      }
      else{
        pin = false;
      }

      PointMass pm = PointMass(position, pin);
      point_masses.push_back(pm);
    }
  }

  // TODO (Part 1.2): Add springs 
  for(int j = 0; j < num_height_points; ++j){
    for(int i = 0; i < num_width_points; ++i){
      // convenient positions
      int index = j * num_height_points + i;
      int left = j * num_height_points + (i-1);
      int above = (j-1) * num_height_points + i;
      int diag_upper_left = (j-1) * num_height_points + (i-1);
      int diag_upper_right = (j-1) * num_height_points + (i+1);
      int two_away_left = j * num_height_points + (i-2);
      int two_away_above = (j-2) * num_height_points + i;

      // structural constraints
      // left
      if(checkValid(i-1, j, num_width_points, num_height_points)){
        Spring sp = Spring(&point_masses[index], &point_masses[left], STRUCTURAL);
        springs.push_back(sp);
      }
      // above
      if(checkValid(i, j-1, num_width_points, num_height_points)){
        Spring sp = Spring(&point_masses[index], &point_masses[above], STRUCTURAL);
        springs.push_back(sp);
      }

      // shearing constraints
      // diag_upper_left
      if(checkValid(i-1, j-1, num_width_points, num_height_points)){
        Spring sp = Spring(&point_masses[index], &point_masses[diag_upper_left], SHEARING);
        springs.push_back(sp);
      }

      // diag_upper_right
      if(checkValid(i+1, j-1, num_width_points, num_height_points)){
        Spring sp = Spring(&point_masses[index], &point_masses[diag_upper_right], SHEARING);
        springs.push_back(sp);
      }

      // bending constraints
      // two_away_left
      if(checkValid(i-2, j, num_width_points, num_height_points)){
        Spring sp = Spring(&point_masses[index], &point_masses[two_away_left], BENDING);
        springs.push_back(sp);
      }

      // two_away_above
      if(checkValid(i, j-2, num_width_points, num_height_points)){
        Spring sp = Spring(&point_masses[index], &point_masses[two_away_above], BENDING);
        springs.push_back(sp);
      }
    }
  }


}

void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * cp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2.1): Compute total force acting on each point mass.

  //Vector3D total_force = Vector3D(0,0,0);
  Vector3D external_force = Vector3D(0,0,0);

  // external force
  for(Vector3D &a: external_accelerations){
    external_force += a;
  }
  external_force *= mass;

  for(PointMass &pm : point_masses){
    pm.forces = external_force;
  }

  // internal force
  for(Spring &s: springs){
    Vector3D direction = s.pm_b->position - s.pm_a->position;
    double length = direction.norm();
    direction.normalize();
    if(cp->enable_structural_constraints && s.spring_type == STRUCTURAL){  
      Vector3D internal_force = cp->ks * (length - s.rest_length) * direction;
      s.pm_a->forces += internal_force;
      s.pm_b->forces += -internal_force;
    }
    
    if(cp->enable_shearing_constraints && s.spring_type == SHEARING){
      Vector3D internal_force = cp->ks * (length - s.rest_length) * direction;
      s.pm_a->forces += internal_force;
      s.pm_b->forces += -internal_force;
    }
   
    if(cp->enable_bending_constraints && s.spring_type == BENDING){
      Vector3D internal_force = 0.2 * cp->ks * (length - s.rest_length) * direction;
      s.pm_a->forces += internal_force;
      s.pm_b->forces += -internal_force;
    }
  }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass positions
    for(PointMass &pm : point_masses){
      if(!pm.pinned){
        Vector3D a = pm.forces / mass;
        Vector3D new_pos = pm.position + (1.0 - cp->damping/100.0) * (pm.position - pm.last_position) + a * delta_t * delta_t;
        pm.last_position = pm.position;
        pm.position = new_pos;
      }
    }


  // This won't do anything until you complete Part 4.
  build_spatial_map();
  for (PointMass &pm : point_masses) {
    self_collide(pm, simulation_steps);
  }

  // This won't do anything until you complete Part 3.
  for (PointMass &pm : point_masses) {
    for (CollisionObject *co : *collision_objects) {
      co->collide(pm);
    }
  }


  // TODO (Part 2.3): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].
  for(Spring &s : springs){
    Vector3D direction = s.pm_b->position - s.pm_a->position;
    double length = direction.norm();
    double distance = length - s.rest_length * 1.1;
    if(distance > 0){
      if(s.pm_a->pinned && !s.pm_b->pinned){
        direction.normalize();
        direction *= distance;
        s.pm_b->position -= direction;
      }

      else if(!s.pm_a->pinned && s.pm_b->pinned){
        direction.normalize();
        direction *= distance;
        s.pm_a->position += direction;
      }

      else if(!s.pm_a->pinned && !s.pm_b->pinned){
        distance /= 2.0;
        direction.normalize();
        direction *= distance;
        s.pm_a->position += direction;
        s.pm_b->position -= direction;
      }
    }
  }
}

void Cloth::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  // TODO (Part 4.2): Build a spatial map out of all of the point masses.
  for(PointMass &pm : point_masses){
    double hash_value = hash_position(pm.position);
    if(!map.count(hash_value)){
      map[hash_value] = new std::vector<PointMass *>();
    }
    map[hash_value]->push_back(&pm);
  }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
  // TODO (Part 4.3): Handle self-collision for a given point mass.
  Vector3D correction = Vector3D(0,0,0);
  float hash_key = hash_position(pm.position);
  double count = 0.0;
  if(map.count(hash_key)){
    vector<PointMass *> *v = map[hash_key];
    for(PointMass *p : *v){
      if(!(p->position == pm.position)){
        double distance = (p->position - pm.position).norm();
        if(distance < 2.0 * thickness){
          double length = 2.0 * thickness - distance;
          Vector3D direction = pm.position - p->position;
          direction.normalize();
          direction *= length;
          correction += direction;
          count += 1.0;
        }
      }
    }
  }
  if(count != 0.0){
    correction = correction / count / simulation_steps;
    pm.position += correction;
  }
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.
  float w = 3.0 * width / num_width_points;
  float h = 3.0 * height / num_height_points;
  float t = max(w, h);

  float x = (pos.x - fmod(pos.x, w)) / w;
  float y = (pos.y - fmod(pos.y, h)) / h;
  float z = (pos.z - fmod(pos.z, t)) / t;

  return (float)(pow(x, 1) + pow(y, 2) + pow(z, 3));
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm, pm + num_width_points, pm + 1));
      triangles.push_back(new Triangle(pm + 1, pm + num_width_points,
                                       pm + num_width_points + 1));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
