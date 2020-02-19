#include "rasterizer.h"
#include <iostream>

using namespace std;
using namespace CGL;
namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                                       size_t width, size_t height,
                                       unsigned int sample_rate) {
  this->psm = psm;
  this->lsm = lsm;
  this->width = width;
  this->height = height;
  this->sample_rate = sample_rate;

  supersample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
  // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
  // NOTE: You are not required to implement proper supersampling for points and lines
  // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
    
  //rgb_framebuffer_target[3 * (y * width + x)    ] = (unsigned char)(c.r * 255);
  //rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(c.g * 255);
  //rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(c.b * 255);
    int i = 0;
    while (i < sample_rate) {
        fill_supersample(x, y, i, c);
        i++;
    }
}

// Optional helper function to add a sample to the supersample_buffer
void RasterizerImp::fill_supersample(size_t x, size_t y, size_t s, Color c) {
  // TODO: Task 2: You may want to implement this function. Hint: our solution uses one line
    supersample_buffer[width * y*sample_rate + x*sample_rate + s] = c;
    //supersample_buffer[width * (y * sample_rate) + (x*sample_rate)] = c;
};

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  fill_pixel(sx, sy, color);
  return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0,
  float x1, float y1,
  Color color) {
  if (x0 > x1) {
    swap(x0, x1); swap(y0, y1);
  }

  float pt[] = { x0,y0 };
  float m = (y1 - y0) / (x1 - x0);
  float dpt[] = { 1,m };
  int steep = abs(m) > 1;
  if (steep) {
    dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
    dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
  }

  while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
    rasterize_point(pt[0], pt[1], color);
    pt[0] += dpt[0]; pt[1] += dpt[1];
  }
}

Vector2D createTangentLines(float x0, float y0, float x1, float y1) {
    Vector2D t;
   
    t.x = x1 - x0;
    t.y = y1 - y0;
    
    return t;
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0,
                                       float x1, float y1,
                                       float x2, float y2,
                                       Color color) {
  // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    Vector2D normCa;
    normCa.x = -1*(y1 - y0);
    normCa.y = x1 - x0;
    Vector2D normCb;
    normCb.x = -1*(y2 - y1);
    normCb.y = x2 - x1;
    Vector2D normCc;
    normCc.x = -1*(y0 - y2);
    normCc.y = x0 - x2;
    
    Vector2D normCCa;
    normCCa.x = -1*(y1 - y2);
    normCCa.y = x1 - x2;
    Vector2D normCCb;
    normCCb.x = -1*(y2 - y0);
    normCCb.y = x2 - x0;
    Vector2D normCCc;
    normCCc.x = -1*(y0 - y1);
    normCCc.y = x0 - x1;
    
    double minX = min(x0, min(x1, x2));
    double minY = min(y0, min(y1, y2));
    double maxX = max(x0, max(x1, x2));
    double maxY = max(y0, max(y1, y2));
    
    for (int i = minX; i < maxX; i++) {
        for (int j = minY; j < maxY; j++) {
            int counter = 0;
            for (int xGrid = 0 ; xGrid < sqrt(sample_rate); xGrid++) {
                for (int yGrid = 0; yGrid < sqrt(sample_rate); yGrid++) {
                    float offSetX;
                    float offSetY;
                    if (sample_rate == 1) {
                        offSetX = 0.5;
                        offSetY = 0.5;
                    } else {
                        offSetX = 1/ sample_rate + xGrid / sqrt(sample_rate);
                        offSetY = 1/ sample_rate + yGrid / sqrt(sample_rate);
                    }
                    Vector2D pointA = createTangentLines(x0, y0, i + offSetX, j + offSetY);
                    Vector2D pointB = createTangentLines(x1, y1, i + offSetX, j + offSetY);
                    Vector2D pointC = createTangentLines(x2, y2, i + offSetX, j + offSetY);
                    
                    float aCROSSbp = normCa.x * pointA.x + pointA.y * normCa.y;
                    float cCROSSap = normCb.x * pointB.x + pointB.y * normCb.y;
                    float bCROSScp = normCc.x * pointC.x + pointC.y * normCc.y;
                    
                    float aCROSSbpCC = normCCa.x * pointC.x + pointC.y * normCCa.y;
                    float cCROSSapCC = normCCb.x * pointA.x + pointA.y * normCCb.y;
                    float bCROSScpCC = normCCc.x * pointB.x + pointB.y * normCCc.y;
                    
                    if ((aCROSSbp >= 0.0f && cCROSSap >= 0.0f && bCROSScp >= 0.0f) || (aCROSSbpCC >= 0.0f && cCROSSapCC >= 0.0f && bCROSScpCC >= 0.0f)) {
                        fill_supersample(i + offSetX, j + offSetY, counter, color);
                    }
                    counter++;
                }
            }
            
        }
    }
            
  // TODO: Task 2: Update to implement super-sampled rasterization
  
}




void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                          float x1, float y1, Color c1,
                                                          float x2, float y2, Color c2)
{
  // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
  // Hint: You can reuse code from rasterize_triangle
    
    double minX = min(x0, min(x1, x2));
    double minY = min(y0, min(y1, y2));
    double maxX = max(x0, max(x1, x2));
    double maxY = max(y0, max(y1, y2));
    
    for (int i = minX; i < maxX; i++) {
        for (int j = minY; j < maxY; j++) {
            int counter = 0;
            for (int xGrid = 0 ; xGrid < sqrt(sample_rate); xGrid++) {
                for (int yGrid = 0; yGrid < sqrt(sample_rate); yGrid++) {
                    double offSetX;
                    double offSetY;
                    if (sample_rate == 1) {
                        offSetX = 0.5;
                        offSetY = 0.5;
                    } else {
                        offSetX = 1/ sample_rate + xGrid / sqrt(sample_rate) + 0.00001;
                        offSetY = 1/ sample_rate + yGrid / sqrt(sample_rate) + 0.00001;
                    }
                    
                    double alpha = (-1 * (i + offSetX - x1) * (y2 - y1) + (j + offSetY - y1) * (x2 - x1)) /
                                    (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
                    
                    double beta = (-1 * (i + offSetX - x2) * (y0 - y2) + (j + offSetY - y2) * (x0 - x2)) /
                                    (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
                    
                    double gamma = 1 - alpha - beta;
                    
                    if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f) {
                        Color t;
                        t.r = alpha * c0.r + beta * c1.r + gamma * c2.r;
                        t.g = alpha * c0.g + beta * c1.g + gamma * c2.g;
                        t.b = alpha * c0.b + beta * c1.b + gamma * c2.b;
                        fill_supersample(i, j, counter, t);
                        
                    }
                    counter++;
                }
            }
            
        }
    }
    

}
double calculateAlpha(int i, double offSetX, int j, double offSetY, float x0, float y0, float x1, float y1, float x2, float y2) {
    
        return (-1 * (i + offSetX - x1) * (y2 - y1) + (j + offSetY - y1) * (x2 - x1)) / (-1 * (x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
}
double calculateBeta(int i, double offSetX, int j, double offSetY, float x0, float y0, float x1, float y1, float x2, float y2) {
        return (-1 * (i + offSetX - x2) * (y0 - y2) + (j + offSetY - y2) * (x0 - x2)) / (-1 * (x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
}

void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
{
    
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    
    SampleParams sampleP;
    sampleP.psm = this->psm;
    sampleP.lsm = this->lsm;
    
    
    
    double minX = min(x0, min(x1, x2));
    double minY = min(y0, min(y1, y2));
    double maxX = max(x0, max(x1, x2));
    double maxY = max(y0, max(y1, y2));
    
    for (int i = minX; i < maxX; i++) {
        for (int j = minY; j < maxY; j++) {
            int counter = 0;
            for (int xGrid = 0 ; xGrid < sqrt(sample_rate); xGrid++) {
                for (int yGrid = 0; yGrid < sqrt(sample_rate); yGrid++) {
                    double offSetX;
                    double offSetY;
                    if (sample_rate == 1) {
                        offSetX = 0.5;
                        offSetY = 0.5;
                    } else {
                        offSetX = 1/ sample_rate + xGrid / sqrt(sample_rate) + 0.00001;
                        offSetY = 1/ sample_rate + yGrid / sqrt(sample_rate) + 0.00001;
                    }
                    
                    double alpha = calculateAlpha(i, offSetX, j, offSetY, x0, y0, x1, y1, x2, y2);
                    double beta = calculateBeta(i, offSetX, j, offSetY, x0, y0, x1, y1, x2, y2);
                    double gamma = 1 - alpha - beta;
                   
                    double alphaDx = calculateAlpha(i + 1, offSetX, j, offSetY, x0, y0, x1, y1, x2, y2);
                    double betaDx = calculateBeta(i + 1, offSetX, j, offSetY, x0, y0, x1, y1, x2, y2);
                    double gammaDx = 1 - alphaDx - betaDx;
                                    
                    double alphaDy = calculateAlpha(i, offSetX, j + 1, offSetY, x0, y0, x1, y1, x2, y2);
                    double betaDy = calculateBeta(i, offSetX, j + 1, offSetY, x0, y0, x1, y1, x2, y2);
                    double gammaDy = 1 - alphaDy - betaDy;
                    
                    if (alpha >= 0.0f && beta >= 0.0f && gamma >= 0.0f) {
                        Vector2D uv;
                        uv.x = alpha*u0 + beta*u1 + gamma*u2;
                        uv.y = alpha*v0 + beta*v1 + gamma*v2;
                        
                        Vector2D uvDx;
                        uvDx.x = alphaDx*u0 + betaDx*u1 + gammaDx*u2;
                        uvDx.y = alphaDx*v0 + betaDx*v1 + gammaDx*v2;
                        
                        Vector2D uvDy;
                        uvDy.x = alphaDy*u0 + betaDy*u1 + gammaDy*u2;
                        uvDy.y = alphaDy*v0 + betaDy*v1 + gammaDy*v2;
                        
                        sampleP.p_uv = uv;
                        sampleP.p_dx_uv = uvDx;
                        sampleP.p_dy_uv = uvDy;
                        
                        Color c = tex.sample(sampleP);
                        
                        fill_supersample(i + offSetX, j + offSetY, counter, c);
                    }
                    counter++;
                }
            }
            
        }
    }
    
    
    
    
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

}

void RasterizerImp::set_sample_rate(unsigned int rate) {
  // TODO: Task 2: You may want to update this function for supersampling support

  this->sample_rate = rate;
    supersample_buffer.resize(width * height * sample_rate, Color::White);

}


void RasterizerImp::set_framebuffer_target( unsigned char* rgb_framebuffer,
                                                size_t width, size_t height )
{
  // TODO: Task 2: You may want to update this function for supersampling support

  this->width = width;
  this->height = height;
    set_sample_rate(sample_rate);
  this->rgb_framebuffer_target = rgb_framebuffer;
  
}


void RasterizerImp::clear_buffers() {
  // TODO: Task 2: You may want to update this function for supersampling support
  // Hint: With supersampling, you have an additional buffer to take care of

  std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    supersample_buffer.clear();
    set_sample_rate(sample_rate);

}


// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
  // TODO: Task 2: You will likely want to update this function for supersampling support
    for (int x = 0; x < width; x++) {
        for (int y = 0; y < height; y++) {
            
            int counter = 0;
            Color t;
            while (counter < sample_rate) {
                t += supersample_buffer[sample_rate*y * width + sample_rate*x + counter];
                //fill_supersample(x, y, counter, t);
                counter++;
            }
            t.r = t.r / sample_rate;
            t.g = t.g / sample_rate;
            t.b = t.b / sample_rate;
            rgb_framebuffer_target[3 * (y * width + x)    ] = (unsigned char)(t.r * 255);
            rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char)(t.g * 255);
            rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char)(t.b * 255);
            
           
        }
    }
  
}

Rasterizer::~Rasterizer() { }


}// CGL
