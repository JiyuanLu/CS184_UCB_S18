#ifndef CGL_IMAGE_H
#define CGL_IMAGE_H

#include "CGL/color.h"
#include "CGL/spectrum.h"

#include <vector>
#include <string.h>

enum Direction {UP = 265, DOWN = 264, LEFT = 263, RIGHT = 262};

namespace CGL {

/**
 * Image buffer which stores color space values with RGBA pixel layout using
 * 32 bit unsigned integers (8-bits per color channel).
 */
struct ImageBuffer {
  /**
   * Default constructor.
   * The default constructor creates a zero-sized image.
   */
  ImageBuffer() : w(0), h(0) {}

  /**
   * Parameterized constructor.
   * Create an image of given size.
   * \param w width of the image
   * \param h height of the image
   */
  ImageBuffer(size_t w, size_t h) : w(w), h(h) { data.resize(w * h); }

  /**
   * Resize the image buffer.
   * \param w new width of the image
   * \param h new height of the image
   */
  void resize(size_t w, size_t h) {
    this->w = w;
    this->h = h;
    data.resize(w * h);
    clear();
  }

  /**
   * Update the color of a given pixel.
   * \param c color value to be set
   * \param x row of the pixel
   * \param y column of the pixel
   */
  void update_pixel(const Color& c, size_t x, size_t y) {
    // assert(0 <= x && x < w);
    // assert(0 <= y && y < h);
    uint32_t p = 0;
    p += ((uint32_t) (clamp(0.f, 1.f, c.a) * 255)) << 24;
    p += ((uint32_t) (clamp(0.f, 1.f, c.b) * 255)) << 16;
    p += ((uint32_t) (clamp(0.f, 1.f, c.g) * 255)) << 8;
    p += ((uint32_t) (clamp(0.f, 1.f, c.r) * 255));
    data[x + y * w] = p;
  }

  /**
   * If the buffer is empty.
   */
  bool is_empty() { return (w == 0 && h == 0); }

  /**
   * Clear image data.
   */
  void clear() { memset(&data[0], 0, w * h * sizeof(uint32_t)); }

  size_t w; ///< width
  size_t h; ///< height
  std::vector<uint32_t> data;  ///< pixel buffer
};

/**
 * High Dynamic Range image buffer which stores linear space spectrum
 * values with 32 bit floating points.
 */
struct HDRImageBuffer {

  /**
   * Default constructor.
   * The default constructor creates a zero-sized image.
   */
  HDRImageBuffer() : w(0), h(0) { }

  /**
   * Parameterized constructor.
   * Create an image of given size.
   * \param w width of the image
   * \param h height of the image
   */
  HDRImageBuffer(size_t w, size_t h) : w(w), h(h) { data.resize(w * h); }

  /**
   * Resize the image buffer.
   * \param w new width of the image
   * \param h new height of the image
   */
  void resize(size_t w, size_t h) {
    this->w = w;
    this->h = h;
    data.resize(w * h);
    clear();
  }

  /**
   * Update the color of a given pixel.
   * \param s new spectrum value to be set
   * \param x row of the pixel
   * \param y column of the pixel
   */
  void update_pixel(const Spectrum& s, size_t x, size_t y) {
    // assert(0 <= x && x < w);
    // assert(0 <= y && y < h);
    data[x + y * w] = s;
  }

  /**
   * Update the color of a given pixel. Blend new pixel color with current
   * pixel data in the buffer according to the blending factor. The result
   * of this would be buffer[i] = s * r + buffer[i] * (1-r);
   * \param s new spectrum value to be set
   * \param x row of the pixel
   * \param y column of the pixel
   * \param r blending factor
   */
  void update_pixel(const Spectrum& s, size_t x, size_t y, float r) {
    // assert(0 <= x && x < w);
    // assert(0 <= y && y < h);
    data[x + y * w] = s * r + (1 - r) * data[x + y * w];
  }

  /**
   * Tonemap and convert to color space image.
   * \param target target color buffer to store output
   * \param gamma gamma value
   * \param level exposure level adjustment
   * \key   key value to map average tone to (higher means brighter)
   * \why   white point (higher means larger dynamic range)
   */
  void tonemap(ImageBuffer& target,
    float gamma, float level, float key, float wht) {

    // compute global log average luminance!
    float avg = 0;
    for (size_t i = 0; i < w * h; ++i) {
      // the small delta value below is used to avoids singularity
      avg += log(0.0000001f + data[i].illum());
    }
    avg = exp(avg / (w * h));


    // apply on pixels
    float one_over_gamma = 1.0f / gamma;
    float exposure = sqrt(pow(2,level));
    for (size_t y = 0; y < h; ++y) {
      for (size_t x = 0; x < w; ++x) {
        Spectrum s = data[x + y * w];
        float l = s.illum();
        s *= key / avg;
        s *= ((l + 1) / (wht * wht)) / (l + 1);
        float r = pow(s.r * exposure, one_over_gamma);
        float g = pow(s.g * exposure, one_over_gamma);
        float b = pow(s.b * exposure, one_over_gamma);
        target.update_pixel(Color(r, g, b, 1.0), x, y);
      }
    }
  }

  /**
   * Convert the given tile of the buffer to color.
   */
  void toColor(ImageBuffer& target, size_t x0, size_t y0, size_t x1, size_t y1) {

    float gamma = 2.2f;
    float level = 1.0f;
    float one_over_gamma = 1.0f / gamma;
    float exposure = sqrt(pow(2,level));
    for (size_t y = y0; y < y1; ++y) {
      for (size_t x = x0; x < x1; ++x) {
        const Spectrum& s = data[x + y * w];
        float r = pow(s.r * exposure, one_over_gamma);
        float g = pow(s.g * exposure, one_over_gamma);
        float b = pow(s.b * exposure, one_over_gamma);
        target.update_pixel(Color(r, g, b, 1.0), x, y);
      }
    }
  }

  /**
   * If the buffer is empty
   */
  bool is_empty() { return (w == 0 && h == 0); }

  /**
   * Clear image buffer.
   */
  void clear() { memset(&data[0], 0, w * h * sizeof(Spectrum)); }

  size_t w; ///< width
  size_t h; ///< height
  std::vector<Spectrum> data; ///< pixel buffer

}; // class HDRImageBuffer

struct LFImageBuffer {

  /**
   * Default constructor.
   * The default constructor creates a zero-sized image.
   */
  LFImageBuffer() : w(0), h(0), subh(0), subw(0) { }

  /**
   * Parameterized constructor.
   * Create an image of given size.
   * \param w width of the image
   * \param h height of the image
   */
  LFImageBuffer(size_t w, size_t h, size_t subw, size_t subh) : w(w), h(h), subh(subh), subw(subw) { 
    
    data.resize(w * h);
    for (auto &grid : data) {
      grid.resize(subh * subw);
    }
    
  }

  /**
   * Resize the image buffer.
   * \param w new width of the image
   * \param h new height of the image
   */
  void resize(size_t w, size_t h, size_t subh, size_t subw) {
    this->w = w;
    this->h = h;
    this->subw = subw;
    this->subh = subh;
    data.resize(w * h);
    for (auto &grid : data) {
      grid.resize(subh * subw);
    }
    clear();
    for (int i = 0 ; i < 10 ; i ++) {
      alphaimage[i].resize(w * h);
      for (auto & grid: alphaimage[i]) {
        grid.resize(subh*subw);
      }
      refocusimage[i].resize(w , h);
    }
    alphamap.resize(w* h);
  }

  /**
   * Update the color of a given pixel.
   * \param s new spectrum value to be set
   * \param x row of the pixel
   * \param y column of the pixel
   */
  void update_pixel(std::vector<Spectrum>& grid, size_t x, size_t y) {
    // assert(0 <= x && x < w);
    // assert(0 <= y && y < h);
    data[x + y * w] = grid;
  }


  /**
   * Convert the given tile of the buffer to color.
   */
  void toColor(ImageBuffer& target, size_t x0, size_t y0, size_t x1, size_t y1) {

    float gamma = 2.2f;
    float level = 1.0f;
    float one_over_gamma = 1.0f / gamma;
    float exposure = sqrt(pow(2,level));
    for (size_t y = y0; y < y1; ++y) {
      for (size_t x = x0; x < x1; ++x) {
        std::vector<Spectrum> grid = data[x + y * w];
        Spectrum s = Spectrum(0, 0, 0);
        /*
        for (auto &spec : grid) {
          s += spec;
        }
        */
       for (int i = 1; i < subh - 1; i++) {
         for (int j = 1; j < subw - 1; j++) {
           s += grid[i * subw + j];
         }
       }

        s = s/float((subw-2) * (subh-2));
//        const Spectrum& s = data[x + y * w];
        float r = pow(s.r * exposure, one_over_gamma);
        float g = pow(s.g * exposure, one_over_gamma);
        float b = pow(s.b * exposure, one_over_gamma);
        target.update_pixel(Color(r, g, b, 1.0), x, y);
      }
    }
  }

  Spectrum getray(double x, double y, double d, int lenx, int leny) {
    double u,v;
    double wstep = 1.0/subw;
    double hstep = 1.0/subh;
    u = lenx*hstep+hstep/2;
    u = u*2*radius-radius;
    v = leny*wstep+wstep/2;
    v = v*2*radius-radius;
//    double x1 = u + (x-u)/(1-d);
//    double y1 = v + (y-v)/(1-d);
//    double x1 = focal/(focal+d)*(u+x)-u;
//    double y1 = focal/(focal+d)*(v+y)-v;
    double x1, y1;
    double C = w/cw*(1/(focal+d)-1/focal);
    x1 = x+ u*C;
    y1 = y+ v*C;
    
    int x0 = floor(x1), y0 = floor(y1);

//    printf("(%f,%f)", u, v);
    Spectrum ret = Spectrum(0,0,0);
    if (0<=x0&&x0<w-1&&y0>=0&&y0<h-1) {
//      printf("!!\n");
      double dx = x1-x0, dy = y1-y0;

      ret = (1-dx)*(1-dy)*data[x0 + y0 * w][lenx*subw + leny]
          + dx*(1-dy)*data[x0+1+y0*w][lenx*subw+leny]
          + (1-dx)*dy*data[x0+y0*w+w][lenx*subw+leny]
          + dx*dy*data[x0+1+y0*w+w][lenx*subw+leny];
//      printf("%f, %f, %f\n", ret.r, ret.g, ret.b);
    }
    return ret;

  }
  void Refocus(ImageBuffer& target, size_t x0, size_t y0, size_t x1, size_t y1, double d, Vector2D pos) {
    float gamma = 2.2f;
    float level = 1.0f;
    float one_over_gamma = 1.0f / gamma;
    float exposure = sqrt(pow(2,level));
    int posx = floor(pos.x*w);
    int posy = floor(pos.y*h);
    int leftx = posx>=3? posx-3:0;
    int lefty = posy>=3? posy-3:0;
    int rx = posx<w-3? posx+3:w;
    int ry = posy<h-3? posy+3:h;
    int a=0;
    for (int i = leftx; i < rx; i ++) {
      for (int j = lefty; j < ry ; j ++) {
//        printf("%d, %d\n", i, j);
        a += alphamap[i+j*w];
//        printf("%d\n", a);
      }
    }
    a = a/((rx-leftx)*(ry-lefty));
//    printf("%d\n", a);
    for (size_t y = y0; y < y1; ++y) {
      for (size_t x = x0; x < x1; ++x) {

        Spectrum s;
        s = refocusimage[a].data[x+y*w];
        float r = pow(s.r * exposure, one_over_gamma);
        float g = pow(s.g * exposure, one_over_gamma);
        float b = pow(s.b * exposure, one_over_gamma);
        
        target.update_pixel(Color(r, g, b, 1.0), x, y);
      }
    }
  }
  void compute_alphamap() {
    static int computed = 0;
    if (computed) {
      return;
    }
    computed = 1;
    printf("compute_alphamap\n");
    double pitch = 2 * radius/subw;
    int x, y, u, v, px, py, a;
    for (a = 0 ; a < 10 ; a ++) {
      printf("%d\n", a);
      double d = -0.5+0.1*a;
      for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
          Spectrum s = Spectrum(0, 0, 0);
          Spectrum t;
          int i, j;
          for (i = 0 ; i < subh; i ++) {
            for (j = 0 ; j < subw; j ++) {
              t = getray(x, y, d, j,i);
              s += t;
              alphaimage[i][x+w*y][j+subw*i] = t;
            }
          }
  //        printf("\n");
          s = s/float(subw * subh);
          refocusimage[a].data[x+y*w] = s;
          
        }
      }
    }
    printf("finish alphaimage\n");
    for (x = 0 ; x < w; x ++) {
      for (y = 0 ; y < h ; y ++) {
        double mins =1000000000000000;
        int mina = -1;
        for (a = 0 ; a < 10 ; a ++) {
          double square_sum = 0;
          double variance = 0;
          Spectrum spa;
          Vector3D rgb;
          for (u = 0 ; u < subw; u ++) {
            for (v = 0 ; v < subh ; v ++) {
              spa = alphaimage[a][x+w*y][u + v*subw];
              rgb = Vector3D(spa.r, spa.g, spa.b);
              square_sum = square_sum + rgb.norm2();
            }
          }
          square_sum = square_sum/(subh*subw);
          spa = refocusimage[a].data[x+y*w];
          rgb = Vector3D(spa.r, spa.g, spa.b);
          variance = square_sum - rgb.norm2();
          if (variance < mins) {
            mins = variance;
            mina = a;
          }
        }
        alphamap[x+y*w] = mina;
      }
    }
    printf("finish alphamap\n");

  }
  void reAperture(ImageBuffer& target, size_t x0, size_t y0, size_t x1, size_t y1, double d, int a) {
    float gamma = 2.2f;
    float level = 1.0f;
    float wstep = 1.0/subw;
    float hstep = 1.0/subh;
    float one_over_gamma = 1.0f / gamma;
    float exposure = sqrt(pow(2,level));
    for (size_t y = y0; y < y1; ++y) {
      for (size_t x = x0; x < x1; ++x) {
        Spectrum s = Spectrum(0, 0, 0);
        int i, j;
        if (a < 10 || a > -10) {
          for (i = 1 - a; i < subh - 1 + a; i ++) {
            for (j = 1 - a ; j < subw - 1 + a; j ++) {
              s += getray(x, y, d, j,i);
            }
          }
          s = s/float((subw-2+a*2) * (subh-2+a*2));
        }
        else {
          //printf("%d\n", a);
          switch(a) {
            case Direction::UP:
              for (i = 4; i < subh; i++) {
                for (j = 2; j < subw - 2; j++) {
                  s+= getray(x, y, d, j, i);
                }
              }
              break;
            case Direction::DOWN:
              for (i = 0; i < subh - 4; i++) {
                for (j = 2; j < subw - 2; j++) {
                  s+= getray(x, y, d, j, i);
                }
              }
              break;
            case Direction::LEFT:
              for (i = 2; i < subh - 2; i++) {
                for (j = 4; j < subw; j++) {
                  s+= getray(x, y, d, j, i);
                }
              }
              break;
            case Direction::RIGHT:
              for (i = 2; i < subh - 2; i++) {
                for (j = 0; j < subw - 4; j++) {
                  s+= getray(x, y, d, j, i);
                }
              }
              break;
          }
          s = s/float((subw-4) * (subh-4));
        }
//        printf("\n");
        float r = pow(s.r * exposure, one_over_gamma);
        float g = pow(s.g * exposure, one_over_gamma);
        float b = pow(s.b * exposure, one_over_gamma);
        target.update_pixel(Color(r, g, b, 1.0), x, y);
      }
    }
  }

  /**
   * If the buffer is empty
   */
  bool is_empty() { return (w == 0 && h == 0); }

  /**
   * Clear image buffer.
   */
  void clear() { 
    for (auto &grid: data) {

      memset(&grid[0], 0, subw*subh*sizeof(Spectrum));

    }
  }

  size_t w; ///< width
  size_t h; ///< height
  size_t subw;
  size_t subh;
  double radius;
  std::vector<std::vector<Spectrum>> data; ///< pixel buffer

  double cw, ch;
  double focal;
  std::vector<std::vector<Spectrum>> alphaimage[10];
  HDRImageBuffer refocusimage[10];
  std::vector<int> alphamap;
}; // class HDRImageBuffer

} // namespace CGL





#endif // CGL_IMAGE_H
