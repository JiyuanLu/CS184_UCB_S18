#include "pathtracer.h"
#include "bsdf.h"
#include "ray.h"

// #include "lenscamera.h"

#include <stack>
#include <random>
#include <algorithm>
#include <sstream>

#include "CGL/CGL.h"
#include "CGL/vector3D.h"
#include "CGL/matrix3x3.h"
#include "CGL/lodepng.h"

#include "GL/glew.h"

#include "static_scene/sphere.h"
#include "static_scene/triangle.h"
#include "static_scene/light.h"

using namespace CGL::StaticScene;

using std::min;
using std::max;

namespace CGL {

PathTracer::PathTracer(size_t ns_aa,
                       size_t max_ray_depth,
                       size_t ns_area_light,
                       size_t ns_diff,
                       size_t ns_glsy,
                       size_t ns_refr,
                       size_t num_threads,
                       size_t samples_per_batch,
                       float max_tolerance,
                       HDRImageBuffer *envmap,
                       bool direct_hemisphere_sample,
                       string filename,
                       double lensRadius,
                       double focalDistance) {
  state = INIT,
      this->ns_aa = ns_aa;
  this->max_ray_depth = max_ray_depth;
  this->ns_area_light = ns_area_light;
  this->ns_diff = ns_diff;
  this->ns_glsy = ns_diff;
  this->ns_refr = ns_refr;
  this->samplesPerBatch = samples_per_batch;
  this->maxTolerance = max_tolerance;
  this->lensRadius = lensRadius;
  this->focalDistance = focalDistance;
  this->direct_hemisphere_sample = direct_hemisphere_sample;
  this->filename = filename;

  if (envmap) {
    this->envLight = new EnvironmentLight(envmap);
  } else {
    this->envLight = NULL;
  }

  bvh = NULL;
  scene = NULL;
  camera = NULL;

  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  show_rays = true;

  imageTileSize = 32;
  numWorkerThreads = num_threads;
  workerThreads.resize(numWorkerThreads);

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;

}

PathTracer::~PathTracer() {

  delete bvh;
  delete gridSampler;
  delete hemisphereSampler;

}

void PathTracer::set_scene(Scene *scene) {

  if (state != INIT) {
    return;
  }

  if (this->scene != nullptr) {
    delete scene;
    delete bvh;
    selectionHistory.pop();
  }

  if (this->envLight != nullptr) {
    scene->lights.push_back(this->envLight);
  }

  this->scene = scene;
  build_accel();

  if (has_valid_configuration()) {
    state = READY;
  }
}

void PathTracer::set_camera(Camera *camera) {
  if (state != INIT) {
    return;
  }
  this->camera = camera;
  if (!this->camera->lensRadius)
    this->camera->lensRadius = lensRadius;
  if (!this->camera->focalDistance)
    this->camera->focalDistance = focalDistance;
  if (has_valid_configuration()) {
    state = READY;
  }
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  if (state != INIT && state != READY) {
    stop();
  }
  sampleBuffer.resize(width, height);
  frameBuffer.resize(width, height);
  cell_tl = Vector2D(0, 0);
  cell_br = Vector2D(width, height);
  render_cell = false;
  sampleCountBuffer.resize(width * height);
  if (has_valid_configuration()) {
    state = READY;
  }
}

bool PathTracer::has_valid_configuration() {
  return scene && camera && gridSampler && hemisphereSampler &&
      (!sampleBuffer.is_empty());
}

void PathTracer::update_screen() {
  switch (state) {
    case INIT:
    case READY:break;
    case VISUALIZE:visualize_accel();
      break;
    case RENDERING:
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      if (render_cell)
        visualize_cell();
      break;
    case DONE:
      //sampleBuffer.tonemap(frameBuffer, tm_gamma, tm_level, tm_key, tm_wht);
      glDrawPixels(frameBuffer.w, frameBuffer.h, GL_RGBA,
                   GL_UNSIGNED_BYTE, &frameBuffer.data[0]);
      if (render_cell)
        visualize_cell();
      break;
  }
}

void PathTracer::stop() {
  switch (state) {
    case INIT:
    case READY:break;
    case VISUALIZE:
      while (selectionHistory.size() > 1) {
        selectionHistory.pop();
      }
      state = READY;
      break;
    case RENDERING:continueRaytracing = false;
    case DONE:
      for (int i = 0; i < numWorkerThreads; i++) {
        workerThreads[i]->join();
        delete workerThreads[i];
      }
      state = READY;
      break;
  }
  render_silent = false;
}

void PathTracer::clear() {
  if (state != READY) return;
  delete bvh;
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  selectionHistory.pop();
  sampleBuffer.resize(0, 0);
  frameBuffer.resize(0, 0);
  state = INIT;
  render_cell = false;
}

void PathTracer::start_visualizing() {
  if (state != READY) {
    return;
  }
  state = VISUALIZE;
}

void PathTracer::start_raytracing() {
  if (state != READY) return;

  // Intersection isect;
  // Ray r = camera->center_ray();
  // if (camera->lens_ind >= 0&& bvh->intersect(r, &isect)) {
  //   camera->focus_at(isect.t);
  // }

  rayLog.clear();
  workQueue.clear();

  state = RENDERING;
  continueRaytracing = true;
  workerDoneCount = 0;

  sampleBuffer.clear();
  if (!render_cell) {
    frameBuffer.clear();
    num_tiles_w = sampleBuffer.w / imageTileSize + 1;
    num_tiles_h = sampleBuffer.h / imageTileSize + 1;
    tilesTotal = num_tiles_w * num_tiles_h;
    tilesDone = 0;
    tile_samples.resize(num_tiles_w * num_tiles_h);
    memset(&tile_samples[0], 0, num_tiles_w * num_tiles_h * sizeof(int));

    // populate the tile work queue
    for (size_t y = 0; y < sampleBuffer.h; y += imageTileSize) {
      for (size_t x = 0; x < sampleBuffer.w; x += imageTileSize) {
        workQueue.put_work(WorkItem(x, y, imageTileSize, imageTileSize));
      }
    }
  } else {
    int w = (cell_br - cell_tl).x;
    int h = (cell_br - cell_tl).y;
    int imTS = imageTileSize / 4;
    num_tiles_w = w / imTS + 1;
    num_tiles_h = h / imTS + 1;
    tilesTotal = num_tiles_w * num_tiles_h;
    tilesDone = 0;
    tile_samples.resize(num_tiles_w * num_tiles_h);
    memset(&tile_samples[0], 0, num_tiles_w * num_tiles_h * sizeof(int));

    // populate the tile work queue
    for (size_t y = cell_tl.y; y < cell_br.y; y += imTS) {
      for (size_t x = cell_tl.x; x < cell_br.x; x += imTS) {
        workQueue.put_work(WorkItem(x, y,
                                    min(imTS, (int) (cell_br.x - x)), min(imTS, (int) (cell_br.y - y))));
      }
    }
  }

  bvh->total_isects = 0;
  bvh->total_rays = 0;
  // launch threads
  fprintf(stdout, "[PathTracer] Rendering... ");
  fflush(stdout);
  for (int i = 0; i < numWorkerThreads; i++) {
    workerThreads[i] = new std::thread(&PathTracer::worker_thread, this);
  }
}

void PathTracer::render_to_file(string filename, size_t x, size_t y, size_t dx, size_t dy) {
  if (x == -1) {
    unique_lock<std::mutex> lk(m_done);
    start_raytracing();
    cv_done.wait(lk, [this] { return state == DONE; });
    lk.unlock();
    save_image(filename);
    fprintf(stdout, "[PathTracer] Job completed.\n");
  } else {
    render_cell = true;
    cell_tl = Vector2D(x, y);
    cell_br = Vector2D(x + dx, y + dy);
    ImageBuffer buffer;
    raytrace_cell(buffer);
    save_image(filename, &buffer);
    fprintf(stdout, "[PathTracer] Cell job completed.\n");
  }
}

void PathTracer::build_accel() {

  // collect primitives //
  fprintf(stdout, "[PathTracer] Collecting primitives... ");
  fflush(stdout);
  timer.start();
  vector<Primitive *> primitives;
  for (SceneObject *obj : scene->objects) {
    const vector<Primitive *> &obj_prims = obj->get_primitives();
    primitives.reserve(primitives.size() + obj_prims.size());
    primitives.insert(primitives.end(), obj_prims.begin(), obj_prims.end());
  }
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // build BVH //
  fprintf(stdout, "[PathTracer] Building BVH from %lu primitives... ", primitives.size());
  fflush(stdout);
  timer.start();
  bvh = new BVHAccel(primitives);
  timer.stop();
  fprintf(stdout, "Done! (%.4f sec)\n", timer.duration());

  // initial visualization //
  selectionHistory.push(bvh->get_root());
}

void PathTracer::visualize_accel() const {

  glPushAttrib(GL_ENABLE_BIT);
  glDisable(GL_LIGHTING);
  glLineWidth(1);
  glEnable(GL_DEPTH_TEST);

  // hardcoded color settings
  Color cnode = Color(.5, .5, .5, .25);
  Color cnode_hl = Color(1., .25, .0, .6);
  Color cnode_hl_child = Color(1., 1., 1., .6);

  Color cprim_hl_left = Color(.6, .6, 1., 1);
  Color cprim_hl_right = Color(.8, .8, 1., 1);
  Color cprim_hl_edges = Color(0., 0., 0., 0.5);

  BVHNode *selected = selectionHistory.top();

  // render solid geometry (with depth offset)
  glPolygonOffset(1.0, 1.0);
  glEnable(GL_POLYGON_OFFSET_FILL);

  if (selected->isLeaf()) {
    bvh->draw(selected, cprim_hl_left);
  } else {
    bvh->draw(selected->l, cprim_hl_left);
    bvh->draw(selected->r, cprim_hl_right);
  }

  glDisable(GL_POLYGON_OFFSET_FILL);

  // draw geometry outline
  bvh->drawOutline(selected, cprim_hl_edges);

  // keep depth buffer check enabled so that mesh occluded bboxes, but
  // disable depth write so that bboxes don't occlude each other.
  glDepthMask(GL_FALSE);

  // create traversal stack
  stack<BVHNode *> tstack;

  // push initial traversal data
  tstack.push(bvh->get_root());

  // draw all BVH bboxes with non-highlighted color
  while (!tstack.empty()) {

    BVHNode *current = tstack.top();
    tstack.pop();

    current->bb.draw(cnode);
    if (current->l) tstack.push(current->l);
    if (current->r) tstack.push(current->r);
  }

  // draw selected node bbox and primitives
  if (selected->l) selected->l->bb.draw(cnode_hl_child);
  if (selected->r) selected->r->bb.draw(cnode_hl_child);

  glLineWidth(3.f);
  selected->bb.draw(cnode_hl);

  // now perform visualization of the rays
  if (show_rays) {
    glLineWidth(1.f);
    glBegin(GL_LINES);

    for (size_t i = 0; i < rayLog.size(); i += 500) {

      const static double VERY_LONG = 10e4;
      double ray_t = VERY_LONG;

      // color rays that are hits yellow
      // and rays this miss all geometry red
      if (rayLog[i].hit_t >= 0.0) {
        ray_t = rayLog[i].hit_t;
        glColor4f(1.f, 1.f, 0.f, 0.1f);
      } else {
        glColor4f(1.f, 0.f, 0.f, 0.1f);
      }

      Vector3D end = rayLog[i].o + ray_t * rayLog[i].d;

      glVertex3f(rayLog[i].o[0], rayLog[i].o[1], rayLog[i].o[2]);
      glVertex3f(end[0], end[1], end[2]);
    }
    glEnd();
  }

  glDepthMask(GL_TRUE);
  glPopAttrib();
}

void PathTracer::visualize_cell() const {
  glPushAttrib(GL_VIEWPORT_BIT);
  glViewport(0, 0, sampleBuffer.w, sampleBuffer.h);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho(0, sampleBuffer.w, sampleBuffer.h, 0, 0, 1);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glTranslatef(0, 0, -1);

  glColor4f(1.0, 0.0, 0.0, 0.8);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_LIGHTING);

  // Draw the Red Rectangle.
  glBegin(GL_LINE_LOOP);
  glVertex2f(cell_tl.x, sampleBuffer.h - cell_br.y);
  glVertex2f(cell_br.x, sampleBuffer.h - cell_br.y);
  glVertex2f(cell_br.x, sampleBuffer.h - cell_tl.y);
  glVertex2f(cell_tl.x, sampleBuffer.h - cell_tl.y);
  glEnd();

  glMatrixMode(GL_PROJECTION);
  glPopMatrix();

  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glPopAttrib();

  glEnable(GL_LIGHTING);
  glEnable(GL_DEPTH_TEST);
}

void PathTracer::key_press(int key) {

  BVHNode *current = selectionHistory.top();
  switch (key) {
    case ']':ns_aa *= 2;
      fprintf(stdout, "[PathTracer] Samples per pixel changed to %lu\n", ns_aa);
      //tm_key = clamp(tm_key + 0.02f, 0.0f, 1.0f);
      break;
    case '[':
      //tm_key = clamp(tm_key - 0.02f, 0.0f, 1.0f);
      ns_aa /= 2;
      if (ns_aa < 1) ns_aa = 1;
      fprintf(stdout, "[PathTracer] Samples per pixel changed to %lu\n", ns_aa);
      break;
    case '=':
    case '+':ns_area_light *= 2;
      fprintf(stdout, "[PathTracer] Area light sample count increased to %zu.\n", ns_area_light);
      break;
    case '-':
    case '_':if (ns_area_light > 1) ns_area_light /= 2;
      fprintf(stdout, "[PathTracer] Area light sample count decreased to %zu.\n", ns_area_light);
      break;
    case '.':
    case '>':max_ray_depth++;
      fprintf(stdout, "[PathTracer] Max ray depth increased to %zu.\n", max_ray_depth);
      break;
    case ',':
    case '<':if (max_ray_depth) max_ray_depth--;
      fprintf(stdout, "[PathTracer] Max ray depth decreased to %zu.\n", max_ray_depth);
      break;
    case ';':
    case ':':focalDistance += .1;
      camera->focalDistance = focalDistance;
      fprintf(stdout, "[PathTracer] Focal distance increased to %f.\n", camera->focalDistance);
      break;
    case '\'':
    case '\"':focalDistance -= .1;
      camera->focalDistance = focalDistance;
      fprintf(stdout, "[PathTracer] Focal distance decreased to %f.\n", camera->focalDistance);
      break;
    case 'k':
    case 'K':
      if (lensRadius == 0)
        lensRadius = .03125f;
      else
        lensRadius *= sqrt(2.);
      camera->lensRadius = lensRadius;
      fprintf(stdout, "[PathTracer] Aperture increased to %f.\n", camera->lensRadius);
      break;
    case 'l':
    case 'L':
      if (lensRadius <= .03125f)
        lensRadius = 0.;
      else
        lensRadius /= sqrt(2.);
      camera->lensRadius = lensRadius;
      fprintf(stdout, "[PathTracer] Aperture decreased to %f.\n", camera->lensRadius);
      break;
    case KEYBOARD_UP:
      if (current != bvh->get_root()) {
        selectionHistory.pop();
      }
      break;
    case KEYBOARD_LEFT:
      if (current->l) {
        selectionHistory.push(current->l);
      }
      break;
    case KEYBOARD_RIGHT:
      if (current->l) {
        selectionHistory.push(current->r);
      }
      break;

    case 'C':render_cell = !render_cell;
      if (render_cell)
        fprintf(stdout, "[PathTracer] Now in cell render mode.\n");
      else
        fprintf(stdout, "[PathTracer] No longer in cell render mode.\n");
      break;

    default:return;
  }
}

Spectrum PathTracer::estimate_direct_lighting_hemisphere(const Ray &r, const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D &hit_p = r.o + r.d * isect.t;
  const Vector3D &w_out = w2o * (-r.d);

  // This is the same number of total samples as estimate_direct_lighting_importance (outside of delta lights).
  // We keep the same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Spectrum L_out = Spectrum(0, 0, 0);

  // TODO (Part 3): Write your sampling loop here
  // COMMENT OUT `normal_shading` IN `est_radiance_global_illumination` BEFORE YOU BEGIN
  Intersection in;
  BVHNode *ro = bvh->get_root();
  //printf("min:%f  %f  %f  max:%f  %f  %f\n", ro->bb.min.x, ro->bb.min.y, ro->bb.min.z, ro->bb.max.x, ro->bb.max.y, ro->bb.max.z);
  for (int i = 0; i < num_samples; i++) {
    Vector3D wi = hemisphereSampler->get_sample();
    float pdf = 1.0 / (2 * PI);

    Vector3D wi_world = o2w * wi;
    Vector3D offset_origin = EPS_D * wi_world + hit_p;
    //printf("origin:%f  %f  %f\n", offset_origin.x, offset_origin.y, offset_origin.z);
    //printf("direct:%f  %f  %f\n", wi_world.x, wi_world.y, wi_world.z);
    Ray out_ray = Ray(offset_origin, wi_world);
    //printf("min_t:%f\n", out_ray.min_t);
    if (bvh->intersect(out_ray, &in, bvh->get_root())) {
      L_out += in.bsdf->get_emission() * cos_theta(wi) * isect.bsdf->f(w_out, wi) / pdf;
      //printf("%f\n", in.bsdf->get_emission().r);
    }
  }

  return L_out / (float) num_samples;
}

Spectrum PathTracer::estimate_direct_lighting_importance(const Ray &r, const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D &hit_p = r.o + r.d * isect.t;
  const Vector3D &w_out = w2o * (-r.d);
  Spectrum L_out = Spectrum(0, 0, 0);

  // TODO (Part 3): Here is where your code for looping over scene lights goes
  // COMMENT OUT `normal_shading` IN `est_radiance_global_illumination` BEFORE YOU BEGIN
  Vector3D wi;
  float distToLight;
  float pdf;
  for (SceneLight *sl : scene->lights) {
    size_t sample_num;
    if (sl->is_delta_light())
      sample_num = 1;
    else
      sample_num = ns_area_light;
    Spectrum temp = Spectrum(0, 0, 0);
    for (int i = 0; i < sample_num; i++) {
      Spectrum radiance = sl->sample_L(hit_p, &wi, &distToLight, &pdf);
      Vector3D w_in = w2o * wi;
      w_in.normalize();
      if (w_in.z < 0)
        continue;
      Ray shadow_ray = Ray(hit_p + wi * EPS_D, wi);
      shadow_ray.max_t = distToLight;
      Intersection is;
      if (bvh->intersect(shadow_ray, &is, bvh->get_root()))
        continue;
      temp += radiance * w_in.z * isect.bsdf->f(w_out, w_in) / pdf;
    }
    L_out += temp / ((float) sample_num);
  }
  return L_out;
}

Spectrum PathTracer::zero_bounce_radiance(const Ray &r, const Intersection &isect) {
  // TODO: Part 4, Task 2
  // Returns the light that results from no bounces of light

  return isect.bsdf->get_emission();
}

Spectrum PathTracer::one_bounce_radiance(const Ray &r, const Intersection &isect) {
  // TODO: Part 4, Task 2
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
  // (you implemented these functions in Part 3)
  if (direct_hemisphere_sample)
    return estimate_direct_lighting_hemisphere(r, isect);
  return estimate_direct_lighting_importance(r, isect);

}

Spectrum PathTracer::at_least_one_bounce_radiance(const Ray &r, const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Spectrum L_out = one_bounce_radiance(r, isect);
  if (max_ray_depth == 0)
    return Spectrum(0, 0, 0);
  else if (max_ray_depth == 1)
    return L_out;
  if (isect.bsdf->is_delta())
    L_out = Spectrum(0, 0, 0);

  // TODO (Part 4.2): Here is where your code for sampling the BSDF,
  // performing Russian roulette step, and returning a recursively
  // traced ray (when applicable) goes
  Vector3D w_in;
  float pdf;
  Spectrum temp = isect.bsdf->sample_f(w_out, &w_in, &pdf);
  float cpdf = 0.7;
  bool ran = coin_flip(cpdf);
  if (ran || r.depth < max_ray_depth) {
    if (r.depth < max_ray_depth)
      cpdf = 1;
    Vector3D wi = o2w * w_in;
    wi.normalize();
    Ray ray = Ray(hit_p + wi * EPS_D, wi);
    ray.depth = r.depth + 1;
    Intersection in;
    if (bvh->intersect(ray, &in, bvh->get_root())) {
      Spectrum L_in = at_least_one_bounce_radiance(ray, in);
      if (isect.bsdf->is_delta())
        L_in += zero_bounce_radiance(ray, in);
      if (pdf > 0)
        L_out += L_in * temp * fabs(w_in.z) / pdf / cpdf;
    }
    else
      if (envLight)
        L_out += envLight->sample_dir(ray) * temp * fabs(w_in.z) / pdf / cpdf;
  }

  return L_out;

}

Spectrum PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Spectrum L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;

  // This line returns a color depending only on the normal vector
  // to the surface at the intersection point.
  // REMOVE IT when you are ready to begin Part 3.

  //return normal_shading(isect.n);

  // TODO (Part 3): Return the direct illumination.
  //L_out += estimate_direct_lighting_hemisphere(r, isect);
  //L_out += estimate_direct_lighting_importance(r, isect);
  // TODO (Part 4): Accumulate the "direct" and "indirect"
  // parts of global illumination into L_out rather than just direct
  L_out = (zero_bounce_radiance(r, isect) + at_least_one_bounce_radiance(r, isect));
  //L_out = at_least_one_bounce_radiance(r, isect) + one_bounce_radiance(r, isect);

  return L_out;
}

Spectrum PathTracer::raytrace_pixel(size_t x, size_t y) {

  // TODO (Part 1.1):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Spectrum.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;            // total samples to evaluate
  Vector2D origin = Vector2D(x, y);    // bottom left corner of the pixel
  if (num_samples == 1) {
    Vector2D sam = gridSampler->get_sample();
    //Ray ray = camera->generate_ray_for_thin_lens((origin.x + 0.5) / (double) sampleBuffer.w, (origin.y + 0.5) / (double) sampleBuffer.h, sam.x, sam.y);
    Ray ray = camera->generate_ray((origin.x + 0.5) / (double) sampleBuffer.w, (origin.y + 0.5) / (double) sampleBuffer.h);
    sampleCountBuffer[y * sampleBuffer.w + x] += 1;
    return est_radiance_global_illumination(ray);
  } else {
    double s1 = 0, s2 = 0;
    Spectrum total = Spectrum(0.0, 0.0, 0.0);
    for (int i = 0; i < num_samples; i++) {
      Vector2D r = gridSampler->get_sample();
      Vector2D sam = gridSampler->get_sample();
      //Ray ray = camera->generate_ray_for_thin_lens((origin.x + r.x) / sampleBuffer.w, (origin.y + r.y) / sampleBuffer.h, sam.x, sam.y);
      Ray ray = camera->generate_ray((origin.x + 0.5) / (double) sampleBuffer.w, (origin.y + 0.5) / (double) sampleBuffer.h);
      Spectrum temp = est_radiance_global_illumination(ray);
      total += temp;
      s1 += temp.illum();
      s2 += temp.illum() * temp.illum();
      if ((i + 1) % samplesPerBatch == 0) {
        double u = s1 / ((double) i + 1);
        double sigma2 = 1.0 / (double) i * (s2 - s1 * s1 / ((double) i + 1));
        double I = 1.96 * sqrt(sigma2) / sqrt((double) i + 1);
        //if (i == num_samples - 1)
        //printf("I:%f    u:%f\n", I, u);
        if (I <= maxTolerance * u) {
          sampleCountBuffer[y * sampleBuffer.w + x] += i + 1;
          //printf("num_samples:%d    samples:%d    samcount:%d\n", num_samples, i+1, sampleCountBuffer[sampleBuffer.w*y+x]);
          return total / (i + 1);
        }
      }
    }
    sampleCountBuffer[y * sampleBuffer.w + x] = num_samples;
    return total / num_samples;
  }
}

void PathTracer::raytrace_tile(int tile_x, int tile_y,
                               int tile_w, int tile_h) {

  size_t w = sampleBuffer.w;
  size_t h = sampleBuffer.h;

  size_t tile_start_x = tile_x;
  size_t tile_start_y = tile_y;

  size_t tile_end_x = std::min(tile_start_x + tile_w, w);
  size_t tile_end_y = std::min(tile_start_y + tile_h, h);

  size_t tile_idx_x = tile_x / imageTileSize;
  size_t tile_idx_y = tile_y / imageTileSize;
  size_t num_samples_tile = tile_samples[tile_idx_x + tile_idx_y * num_tiles_w];

  for (size_t y = tile_start_y; y < tile_end_y; y++) {
    if (!continueRaytracing) return;
    for (size_t x = tile_start_x; x < tile_end_x; x++) {
      Spectrum s = raytrace_pixel(x, y);
      sampleBuffer.update_pixel(s, x, y);
    }
  }

  tile_samples[tile_idx_x + tile_idx_y * num_tiles_w] += 1;
  sampleBuffer.toColor(frameBuffer, tile_start_x, tile_start_y, tile_end_x, tile_end_y);
}

void PathTracer::raytrace_cell(ImageBuffer &buffer) {
  size_t tile_start_x = cell_tl.x;
  size_t tile_start_y = cell_tl.y;

  size_t tile_end_x = cell_br.x;
  size_t tile_end_y = cell_br.y;

  size_t w = tile_end_x - tile_start_x;
  size_t h = tile_end_y - tile_start_y;
  HDRImageBuffer sb(w, h);
  buffer.resize(w, h);

  stop();
  render_cell = true;
  {
    unique_lock<std::mutex> lk(m_done);
    start_raytracing();
    cv_done.wait(lk, [this] { return state == DONE; });
    lk.unlock();
  }

  for (size_t y = tile_start_y; y < tile_end_y; y++) {
    for (size_t x = tile_start_x; x < tile_end_x; x++) {
      buffer.data[w * (y - tile_start_y) + (x - tile_start_x)] = frameBuffer.data[x + y * sampleBuffer.w];
    }
  }
}

void PathTracer::worker_thread() {

  Timer timer;
  timer.start();

  WorkItem work;
  while (continueRaytracing && workQueue.try_get_work(&work)) {
    raytrace_tile(work.tile_x, work.tile_y, work.tile_w, work.tile_h);
    {
      lock_guard<std::mutex> lk(m_done);
      ++tilesDone;
      if (!render_silent)
        cout << "\r[PathTracer] Rendering... " << int((double) tilesDone / tilesTotal * 100) << '%';
      cout.flush();
    }
  }

  workerDoneCount++;
  if (!continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    if (!render_silent) fprintf(stdout, "\n[PathTracer] Rendering canceled!\n");
    state = READY;
  }

  if (continueRaytracing && workerDoneCount == numWorkerThreads) {
    timer.stop();
    if (!render_silent) fprintf(stdout, "\r[PathTracer] Rendering... 100%%! (%.4fs)\n", timer.duration());
    if (!render_silent) fprintf(stdout, "[PathTracer] BVH traced %llu rays.\n", bvh->total_rays);
    if (!render_silent)
      fprintf(stdout, "[PathTracer] Averaged %f intersection tests per ray.\n",
              (((double) bvh->total_isects) / bvh->total_rays));

    lock_guard<std::mutex> lk(m_done);
    state = DONE;
    cv_done.notify_one();
  }
}

void PathTracer::save_image(string filename, ImageBuffer *buffer) {

  if (state != DONE) return;

  if (!buffer)
    buffer = &frameBuffer;

  if (filename == "") {
    time_t rawtime;
    time(&rawtime);

    time_t t = time(nullptr);
    tm *lt = localtime(&t);
    stringstream ss;
    ss << this->filename << "_screenshot_" << lt->tm_mon + 1 << "-" << lt->tm_mday << "_"
       << lt->tm_hour << "-" << lt->tm_min << "-" << lt->tm_sec << ".png";
    filename = ss.str();
  }

  uint32_t *frame = &buffer->data[0];
  size_t w = buffer->w;
  size_t h = buffer->h;
  uint32_t *frame_out = new uint32_t[w * h];
  for (size_t i = 0; i < h; ++i) {
    memcpy(frame_out + i * w, frame + (h - i - 1) * w, 4 * w);
  }

  fprintf(stderr, "[PathTracer] Saving to file: %s... ", filename.c_str());
  lodepng::encode(filename, (unsigned char *) frame_out, w, h);
  fprintf(stderr, "Done!\n");

  save_sampling_rate_image(filename);
}

void PathTracer::save_sampling_rate_image(string filename) {
  size_t w = frameBuffer.w;
  size_t h = frameBuffer.h;
  ImageBuffer outputBuffer(w, h);

  for (int x = 0; x < w; x++) {
    for (int y = 0; y < h; y++) {
      float samplingRate = sampleCountBuffer[y * w + x] * 1.0f / ns_aa;

      Color c;
      if (samplingRate <= 0.5) {
        float r = (0.5 - samplingRate) / 0.5;
        c = Color(0.0f, 0.0f, 1.0f) * r + Color(0.0f, 1.0f, 0.0f) * (1.0 - r);
      } else {
        float r = (1.0 - samplingRate) / 0.5;
        c = Color(0.0f, 1.0f, 0.0f) * r + Color(1.0f, 0.0f, 0.0f) * (1.0 - r);
      }
      outputBuffer.update_pixel(c, x, h - 1 - y);
    }
  }

  lodepng::encode(filename.substr(0, filename.size() - 4) + "_rate.png",
                  (unsigned char *) (outputBuffer.data.data()), w, h);
}

}  // namespace CGL
