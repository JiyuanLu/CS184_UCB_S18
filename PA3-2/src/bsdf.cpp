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

  return reflectance / PI;
}

Spectrum DiffuseBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO (Part 3.1): 
  // This function takes in only wo and provides pointers for wi and pdf,
  // which should be assigned by this function.
  // After sampling a value for wi, it returns the evaluation of the BSDF
  // at (wo, *wi).
  *wi = sampler.get_sample(pdf);
  return f(wo, *wi);
}


// Mirror BSDF //

Spectrum MirrorBSDF::f(const Vector3D& wo, const Vector3D& wi) {
  return Spectrum();
}

Spectrum MirrorBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
  // TODO:
  // Implement MirrorBSDF
  *pdf = 1.0;
  reflect(wo, wi);
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
  Vector3D n(0,0,1);
  double tan_theta, numerator, denominator;
  tan_theta = sin_theta(h) / cos_theta(h);
  numerator = exp(-pow(tan_theta, 2) / pow(alpha, 2));
  denominator = PI * pow(alpha, 2) * pow(cos_theta(h), 4);
  return numerator / denominator;
}

Spectrum MicrofacetBSDF::F(const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Compute Fresnel term for reflection on dielectric-conductor interface.
    // You will need both eta and etaK, both of which are Spectrum.
    Spectrum Rs, Rp, eta_k;
    double c = cos_theta(wi);
    double c_2 = pow(c, 2);
    eta_k = eta * eta + k * k;
    Rs = (eta_k - 2 * eta * c + c_2) / (eta_k + 2 * eta * c + c_2);
    Rp = (eta_k * c_2 - 2 * eta * c + 1) / (eta_k * c_2 + 2 * eta * c + 1);
    return (Rs + Rp) / 2;
}

Spectrum MicrofacetBSDF::f(const Vector3D& wo, const Vector3D& wi) {
    // TODO: proj3-2, part 2
    // Implement microfacet model here.
    Vector3D n(0,0,1);
    if((dot(n, wo) > 0) && (dot(n, wi) > 0)){
      Vector3D h = wo + wi;
      h.normalize();
      return F(wi) * G(wo, wi) * D(h) / 4.0 / dot(n, wo) / dot(n, wi);
    }
    else
      return Spectrum();
}

Spectrum MicrofacetBSDF::sample_f(const Vector3D& wo, Vector3D* wi, float* pdf) {
    // TODO: proj3-2, part 2
    // *Importance* sample Beckmann normal distribution function (NDF) here.
    // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
    //       and return the sampled BRDF value.

  //*wi = cosineHemisphereSampler.get_sample(pdf); //placeholder
  //return MicrofacetBSDF::f(wo, *wi);

  
  Vector2D random_pair = sampler.get_sample();

  float theta_h = atan(sqrt(-pow(alpha, 2) * log(1.0 - random_pair.x)));
  float phi_h = 2.0 * PI * random_pair.y;
  float h_x = cos(phi_h) * sin(theta_h);
  float h_y = sin(phi_h) * sin(theta_h);
  float h_z = cos(theta_h);
  Vector3D h(h_x, h_y, h_z);
  *wi = -wo + 2.0 * dot(wo, h) * h;

  float pdf_theta_h = 2.0 * sin(theta_h) * exp(- pow(tan(theta_h), 2) / pow(alpha, 2)) / pow(alpha, 2) / pow(cos(theta_h), 3);
  float pdf_phi_h = 1.0 / 2.0 / PI;
  float pdf_sample_h_w = pdf_theta_h * pdf_phi_h / sin(theta_h);
  float pdf_sample_wi_w = pdf_sample_h_w / 4.0 / dot(*wi, h);
  *pdf = pdf_sample_wi_w;
  return MicrofacetBSDF::f(wo, *wi);
  
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
  
  // total internal reflection
  if(!refract(wo,wi,ior)){
    *pdf = 1.0;
    reflect(wo, wi);
    return reflectance / abs_cos_theta(*wi);
  }

  else{
    float R_o = pow(((1.0-ior)/ (1.0+ior)), 2);
    float R = R_o + (1.0 - R_o) * pow((1.0 - abs_cos_theta(wo)), 5);
    if(coin_flip(R)){
      *pdf = R;
      reflect(wo, wi);
      return R * reflectance / abs_cos_theta(*wi);
    }
    else{
      refract(wo,wi,ior);
      *pdf = 1.0 - R;
      float eta;

      if(wo.z >= 0)
        eta = 1.0 / ior;
      else
        eta = ior;

      return (1.0 - R) * transmittance / abs_cos_theta(*wi) / pow(eta, 2);
    }
  }
}

void BSDF::reflect(const Vector3D& wo, Vector3D* wi) {
  // TODO: 3-2 Part 1 Task 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
  *wi = Vector3D(-wo[0], -wo[1], wo[2]);
}

bool BSDF::refract(const Vector3D& wo, Vector3D* wi, float ior) {
  // TODO: 3-2 Part 1 Task 3
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When wo.z is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  float eta, inside;
  float wi_x, wi_y, wi_z;
  //enter
  if(wo.z >= 0){
    eta = 1.0 / ior;
    inside = 1.0 - pow(eta, 2) * (1.0 - pow(wo.z, 2));
    if(inside < 0)
      return false;
    wi_z = - sqrt(inside);
  }

  //exit
  else{
    eta = ior;
    inside = 1.0 - pow(eta, 2) * (1.0 - pow(wo.z, 2));
    if(inside < 0)
      return false;
    wi_z = sqrt(inside);
  }

  if ((1.0 - pow(eta, 2) * (1.0 - pow(wo.z, 2))) < 0)
    return false;

  wi_x = -eta * wo.x;
  wi_y = -eta * wo.y;
  *wi = Vector3D(wi_x, wi_y, wi_z);
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
