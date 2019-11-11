#include "texture.h"
#include "CGL/color.h"

namespace CGL {

Color Texture::sample(const SampleParams &sp) {
  // Part 5: Fill this in.
  if (sp.lsm == L_ZERO){
    if (sp.psm == P_NEAREST){
      return sample_nearest(sp.p_uv,0);
    }
    else{
      return sample_bilinear(sp.p_uv,0);
    }
  }
  else if (sp.lsm == L_NEAREST){
    if(sp.psm == P_NEAREST){
      return sample_nearest(sp.p_uv,get_level(sp));
    }
    else{
      return sample_bilinear(sp.p_uv,0);
    }
  }
  else{
      return sample_trilinear(sp.p_uv, get_level(sp));
  }
}

float Texture::get_level(const SampleParams &sp) {
  // Part 6: Fill this in.
  Vector2D dx = sp.p_dx_uv - sp.p_uv;
  Vector2D dy = sp.p_dy_uv - sp.p_uv;
  dx.x *= this -> width;
  dx.y *= this -> height;
  dy.x *= this -> width;
  dy.y *= this -> height;
  float length = max(sqrt(pow(dx.x,2) + pow(dx.y,2)), sqrt(pow(dy.x,2)) + pow(dy.y,2));
  float level = log2(length);
  if(level > (this -> mipmap.size() -1)){
    level = this -> mipmap.size() - 1;
  }
  else if(level < 0){
    level = 0;
  }
  return level;
}

Color Texture::sample_nearest(Vector2D uv, float level) {
  // Part 5: Fill this in.
  // default level = 0
  level = (int)round(level);
  int w = mipmap[level].width;
  int h = mipmap[level].height;

  int X = uv.x * w;
  int Y = uv.y * h;
  if(X>=0 && X<w && Y>=0 && Y<h){
    int index = 4 * (Y * w + X);
    float r = mipmap[level].texels[index] / 255.0;
    float g = mipmap[level].texels[index+1] / 255.0;
    float b = mipmap[level].texels[index+2] / 255.0;
    return Color(r,g,b);
  }
}

Color lerp(float x, Color c0, Color c1){
  return c0 + x * (c1 - c0);
}

Color Texture::sample_bilinear(Vector2D uv, float level) {
  // Part 5: Fill this in.
  // default level = 0
  level = (int)round(level);
  int w = mipmap[level].width;
  int h = mipmap[level].height;

  int x0 = floor(uv.x * w);
  int y0 = floor(uv.y * h);
  int x1 = ceil(uv.x * w);
  int y1 = ceil(uv.y * h);

  if(x0 >= 0 && x0 < w+1 && x1 >=0 && x1 < w+1 && y0 >= 0 && y0 < h+1 && y1 >=0 && y1 < h+1){
    int texel_index_00 = 4 * (y0 * w + x0);
    int texel_index_01 = 4 * (y0 * w + x1);
    int texel_index_10 = 4 * (y1 * w + x0);
    int texel_index_11 = 4 * (y1 * w + x1);

    Color p00 = Color(mipmap[level].texels[texel_index_00] / 255.0,
                      mipmap[level].texels[texel_index_00 + 1] / 255.0,
                      mipmap[level].texels[texel_index_00 + 2] / 255.0);
    Color p01 = Color(mipmap[level].texels[texel_index_01] / 255.0,
                      mipmap[level].texels[texel_index_01 + 1] / 255.0,
                      mipmap[level].texels[texel_index_01 + 2] / 255.0);
    Color p10 = Color(mipmap[level].texels[texel_index_10] / 255.0,
                      mipmap[level].texels[texel_index_10 + 1] / 255.0,
                      mipmap[level].texels[texel_index_10 + 2] / 255.0);
    Color p11 = Color(mipmap[level].texels[texel_index_11] / 255.0,
                      mipmap[level].texels[texel_index_11 + 1] / 255.0,
                      mipmap[level].texels[texel_index_11 + 2] / 255.0);
    Color u0 = lerp(uv.x - floor(uv.x), p00, p01);
    Color u1 = lerp(uv.x - floor(uv.x), p10, p11);
    return lerp(uv.y - floor(uv.y), u0, u1);
   }
}

//Color Texture::sample_trilinear(Vector2D uv, Vector2D du, Vector2D dv) {
  // Part 6: Fill this in.
  //return Color();
//}

Color Texture::sample_trilinear(Vector2D uv, float level){
  Color down = sample_bilinear(uv,floor(level));
  Color up = sample_bilinear(uv,ceil(level));
  return lerp((level - floor(level)), down, up);
}




/****************************************************************************/



inline void uint8_to_float(float dst[4], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
  dst[3] = src_uint8[3] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[4]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
  dst_uint8[3] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[3])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    //assert (width > 0);
    height = max(1, height / 2);
    //assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(4 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 4; // 32 bit RGBA
    int currLevelPitch = currLevel.width * 4; // 32 bit RGBA

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[4];
    float input[4];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      //assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = result[3] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 4 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
          result[3] += wWeight[ii] * input[3];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (4 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      //assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = result[3] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
          result[3] += hWeight[jj] * input[3];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = result[3] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        4 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
              result[3] += weight * input[3];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 4 * i, result);
        }
      }
    }
  }
}

}
