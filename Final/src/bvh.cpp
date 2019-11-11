#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL { namespace StaticScene {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  root = construct_bvh(_primitives, max_leaf_size);

}

BVHAccel::~BVHAccel() {
  if (root) delete root;
}

BBox BVHAccel::get_bbox() const {
  return root->bb;
}

void BVHAccel::draw(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->draw(c);
  } else {
    draw(node->l, c);
    draw(node->r, c);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color& c) const {
  if (node->isLeaf()) {
    for (Primitive *p : *(node->prims))
      p->drawOutline(c);
  } else {
    drawOutline(node->l, c);
    drawOutline(node->r, c);
  }
}

BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox centroid_box, bbox;

  for (Primitive *p : prims) {
    BBox bb = p->get_bbox();
    bbox.expand(bb);
    Vector3D c = bb.centroid();
    centroid_box.expand(c);
  }

  BVHNode *node = new BVHNode(bbox);

  

  if (prims.size()<= max_leaf_size) {
    node->prims = new vector<Primitive *>(prims);
    return node;
  }
  std::vector<Primitive*> lp, rp;

  int a, i;
  double t = -1;
  Vector3D center = centroid_box.centroid();

  for (i = 0 ; i < 3; i ++) {
    if(centroid_box.extent[i]>t) {
      t = centroid_box.extent[i];
      a = i;

    }
  }
  for (Primitive *p : prims) {
    BBox pb = p->get_bbox();
    Vector3D pc = pb.centroid();
    if (pc[a]<center[a]) {
      lp.push_back(p);
    }
    else {
      rp.push_back(p);
    }
  }


  node->l = construct_bvh(lp, max_leaf_size);
  node->r = construct_bvh(rp, max_leaf_size);

  return node;

}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.
  bool hl, hr;
  double t0 = ray.min_t, t1 = ray.max_t;
  if (node->bb.intersect(ray, t0, t1)) {
    if (node->isLeaf()) {
      for (Primitive *p : *(node->prims)) {
        total_isects++;
        if (p->intersect(ray)) 
        return true;
      }
      
    } 
    hl = intersect(ray, node->l);
    hr = intersect(ray, node->r);
    if (hl || hr) 
    return true;
  }
  return false;

}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  
  bool hl, hr, a, b;
  bool hit = false;
  double t0 = ray.min_t, t1 = ray.max_t;
  //printf("%f, %f, %f, %f, %f, %f\n", node->bb.min.x, node->bb.min.y
   // , node->bb.min.z, node->bb.max.x, node->bb.max.y, node->bb.max.z);
  if (node->bb.intersect(ray, t0, t1)) {
    if (node->isLeaf()) {
      for (Primitive *p : *(node->prims)) {
//        printf("ray dir:%f %f %f, time: %f\n", ray.d.x, ray.d.y, ray.d.z, ray.max_t);
        total_isects++;
        if (p->intersect(ray, i)) {
          hit = true;
        }
      }
      return hit;
    }
    
    hl = intersect(ray, i, node->l);
    hr = intersect(ray, i, node->r);
    hit = hl || hr;
    return hit;
  }

  return hit;
  /*
  
  bool hit = false;
  for (Primitive *p : *(root->prims)) {
    total_isects++;
    if (p->intersect(ray, i))
      hit = true;
  }
  return hit;
  */
}

}  // namespace StaticScene
}  // namespace CGL
