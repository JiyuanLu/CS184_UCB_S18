#include "bsdf.h"

#include <iostream>
#include <algorithm>
#include <utility>

using std::min;
using std::max;
using std::swap;

namespace CGL {

void make_coord_space(Matrix3x3& o2w, const Vector3D& n) {

    Vector3D z = Vector3D(n.x, n.y, n.z);
    Vector3D h = z;
    if (fabs(h.x) <= fabs(h.y) && fabs(h.x) <= fabs(h.z)) h.x = 1.0;
    else if (fabs(h.y) <= fabs(h.x) && fabs(h.y) <= fabs(h.z)) h.y = 1.0;
    else h.z = 1.0;

    z.normalize();
    Vector3D y = cross(h, z);
    y.normalize();
    Vector3D x = cross(z, y);
    x.normalize();

    o2w[0] = x;
    o2w[1] = y;
    o2w[2] = z;
}

// Diffuse BSDF //

Spectrum DiffuseBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  // TODO (Part 3.1): 
  // This function takes in both wo and wi and returns the evaluation of
  // the BSDF for those two directions.

  return reflectance*(1.0/PI);
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (Part 3.1): 
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).

  *wi = sampler.get_sample(pdf);

  return reflectance*(1.0/PI);
}


// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement MirrorBSDF
  reflect(wo, wi);
  *pdf = 1;

  return reflectance / abs_cos_theta(*wi);
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D& wo, const Vector3D& wi) {
    return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D& h) {
    // TODO: proj3-2, part 2
    // Compute Beckmann normal distribution function (NDF) here.
    // You will need the roughness alpha.
    double c = cos_theta(h);
    double s = sin_theta(h);
    double t = s/c;
    return exp(-t*t/alpha/alpha)/PI/alpha/alpha/pow(c, 4);
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and etaK, both of which are Spectrum.
    Spectrum rs, rp;
    Spectrum t = eta*eta + k*k;
    double c = fabs(wi.z);
    rs = (t-2*eta*c+ c*c)/(t+2*eta*c + c*c);
    rp = (t*c*c-2*eta*c+1)/(t*c*c+2*eta*c+1);
    return (rs+rp)*0.5;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Implement microfacet model here.
    if (wo.z<0|| wi.z < 0) {
      return Spectrum();
    }
    Vector3D h = wo+wi;
    Vector3D n(0,0,1);
    h.normalize();
    Spectrum res = F(wi)*G(wo, wi)*D(h)/4/wo.z/wi.z;
    return res;
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.
  
    Vector2D r = sampler.get_sample();
    double th2 = -alpha*alpha*log(1-r.x);
    double c = sqrt(1/(1+th2));
    double s = sqrt(1-1/(1+th2));
    Vector3D h;
    h.z = c;
    h.x = cos(2.*PI*r.y)*s;
    h.y = -sin(2.*PI*r.y)*s;
    h.normalize();
//    Vector3D t = Vector3D(wo.x, wo.z, fabs(wo.z));
    *wi = 2*dot(h, wo)*h-wo;
    if (wi->z <0){
      *pdf = 1;
      return Spectrum();
    }
    wi->normalize();
    double ptheta = 2*s/alpha/alpha/pow(c,3)*exp(-th2/alpha/alpha);
    double ph = ptheta/2./PI/s;
    *pdf = ph/4/dot(*wi, h);
    if (*pdf < 0) {
      *pdf = 1;
      return Spectrum();
    }
    return MicrofacetBSDF::f(wo, *wi);
    
  
//    *wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
//    return MicrofacetBSDF::f(wo, *wi);
    
}

// Refraction BSDF //

Spectrum RefractionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum RefractionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {

  // TODO:
  // Implement RefractionBSDF
  return Spectrum();
}

// Glass BSDF //

Spectrum GlassBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum GlassBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO: 3-2 Part 1 Task 4
  // Compute Fresnel coefficient and either reflect or refract based on it.
  if (!refract(wo, wi, ior)){    
    reflect(wo, wi);
    *pdf = 1.0;
    return reflectance / abs_cos_theta(*wi);
  }
  double eta;
  if (wo.z>0)
    eta = 1/ior;
  else 
    eta = ior;
  double R0 = (1-eta)/(1+eta);
  R0 *= R0;
  double R = R0 + (1-R0)*pow(1-fabs(wi->z), 5);
  if (coin_flip(R)) {
    reflect(wo, wi);
    *pdf = R;
    return R * reflectance/ abs_cos_theta(*wi);
  }
  else {
    refract(wo, wi, ior);
    *pdf = 1-R;
    return (1-R) * transmittance / abs_cos_theta(*wi) / eta*eta;
  }
}


void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO: 3-2 Part 1 Task 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  Vector3D normal(0,0,1);
  wi->x = -wo.x;
  wi->y = -wo.y;
  wi->z = wo.z;
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO: 3-2 Part 1 Task 3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When wo.z is positive, then wo corresponds to a
  // ray entering the surface through vacuum.
  double eta;
  if (wo.z>0)
    eta = 1/ior;
  else {
    eta = ior;
  }
  double d = 1-eta*eta*(1-wo.z*wo.z);
  if (d<0)
    return false;
  wi->x = -eta*wo.x;
  wi->y = -eta* wo.y;
  if (wo.z>0)
    wi->z = -sqrt(d);
  else wi->z = sqrt(d);
  return true;
}

// Emission BSDF //

Spectrum EmissionBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum EmissionBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  *pdf = 1.0 / PI;
  *wi  = sampler.get_sample(pdf);
  return Spectrum();
}

} // namespace CGL
