#include "bvh.h"

#include "CGL/CGL.h"
#include "static_scene/triangle.h"

#include <iostream>
#include <stack>
#include <algorithm>

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


int largestDimension(Vector3D v) {
  if (v.x >= v.y && v.x >= v.z) {
    return 0;
  }
  else if (v.y >= v.x && v.y >= v.z) {
    return 1;
  }
  // else {
  else
    return 2;
}


BVHNode *BVHAccel::construct_bvh(const std::vector<Primitive*>& prims, size_t max_leaf_size) {
  
  // Part 2, Task 1:
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
  
  if (prims.size() <= max_leaf_size) {
    node->prims = new vector<Primitive *>(prims);
    return node;
  }

  else {
    int axis = largestDimension(centroid_box.extent);

    vector<Primitive *> prims1 = vector<Primitive *>();
    vector<Primitive *> prims2 = vector<Primitive *>();
    
    bool p1 = false;
    bool p2 = false;
    for (Primitive *p : prims) {
      if (p->get_bbox().centroid()[axis] <= centroid_box.centroid()[axis]) {
        prims1.push_back(p);
        p1 = true;
      }
      else {
        prims2.push_back(p);
        p2 = true;
      }
    }

    if (!p1 || !p2) {
      int half_prims_size = prims.size() / 2;

      vector<Primitive *> left = vector<Primitive *>();
      vector<Primitive *> right = vector<Primitive *>();

      if (!p1) {
        for (int i = 0; i < half_prims_size; i++) {
          left.push_back(prims2[i]);
        }

        for (int i = half_prims_size; i < prims.size(); i++) {
          right.push_back(prims2[i]);
        }
        node->l = construct_bvh(left, max_leaf_size); 
        node->r = construct_bvh(right, max_leaf_size);
      }
      else {
        for (int i = 0; i < half_prims_size; i++) {
          left.push_back(prims1[i]);
        }
        for (int i = half_prims_size; i < prims.size(); i++) {
          right.push_back(prims1[i]);
        }
        node->l = construct_bvh(left, max_leaf_size);    
        node->r = construct_bvh(right, max_leaf_size);
      }
    }
    else {
      node->l = construct_bvh(prims1, max_leaf_size);
      node->r = construct_bvh(prims2, max_leaf_size);
    }
    return node;
  }
}


bool BVHAccel::intersect(const Ray& ray, BVHNode *node) const {
  // Part 2, task 3: replace this.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

  double t0, t1;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  else {
    if (t1 < ray.min_t || t0 > ray.max_t) {
      return false;
    }
    else {
      if (node->isLeaf()) {
        bool hit = false;
        for (Primitive *p : *(node->prims)) {
          total_isects++;
          if (p->intersect(ray)) {
            hit = true;
          }
        }
        return hit;
      }
      else {
        bool left = false;
        bool right = false;
        left = left || intersect(ray, node->l);
        right = right || intersect(ray, node->r);
        return left || right;
      }
    }
  }
}

bool BVHAccel::intersect(const Ray& ray, Intersection* i, BVHNode *node) const {
  // Part 2, task 3: replace this

  double t0, t1;
  if (!node->bb.intersect(ray, t0, t1)) {
    return false;
  }
  else {
    if (t1 < ray.min_t || t0 > ray.max_t) {
      return false;
    }

    else {
      if (node->isLeaf()) {
        bool hit = false;
        for (Primitive *p : *(node->prims)) {
          total_isects++;
          if (p->intersect(ray, i)) {
            hit = true;
          }
        }
        return hit;
      }
      else {
        bool left = false;
        bool right = false;
        left = left || intersect(ray, i, node->l);
        right = right || intersect(ray, i, node->r);
        return left || right;
      } 
    }
  }
}

}  // namespace StaticScene
}  // namespace CGL
