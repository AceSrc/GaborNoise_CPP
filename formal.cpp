/* noise.cpp
 *
 * Example code for the article
 *
 * Procedural Noise using Sparse Gabor Convolution.
 * Ares Lagae, Sylvain Lefebvre, George Drettakis and Philip Dutre.
 * ACM Transactions on Graphics (Proceedings of ACM SIGGRAPH 2009) 28(3), 2009.
 *
 * If you use this code, please reference the article.
 *
 * This code produces an image 'noise.ppm' of anisotropic noise for the
 * parameters in the 'main' function. This image should be similar to the image
 * 'noise_reference.ppm'.
 *
 * This code allows to generate anisotropic/isotropic and periodic/nonperiodic
 * noise. Please refer to the comments in the code.
 *
 * This code is designed for brevity and clarity. It is significantly slower
 * than the GPU version.
 *
 * Please report bugs to the authors.
 *
 * Revision history
 *   2009-05-06  initial version
 */

/* Copyright (c) 2009 by Ares Lagae, Sylvain Lefebvre,
 * George Drettakis and Philip Dutre
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

// -----------------------------------------------------------------------------

#include <climits>
#include <cmath>
#include <iostream>

#ifndef M_PI
#  define M_PI 3.14159265358979323846
#endif

class pseudo_random_number_generator
{
public:
  void seed(unsigned s) { x_ = s; }
  unsigned operator()() { x_ *= 3039177861u; return x_; }
  float uniform_0_1() { return float(operator()()) / float(UINT_MAX); }
  float uniform(float min, float max)
    { return min + (uniform_0_1() * (max - min)); }
  unsigned poisson(float mean)
  {
    float g_ = std::exp(-mean);
    unsigned em = 0;
    double t = uniform_0_1();
    while (t > g_) {
      ++em;
      t *= uniform_0_1();
    }
    return em;
  }
private:
  unsigned x_;
};

float gabor(float K, float a, float F_0, float omega_0, float x, float y)
{
  float gaussian_envelop = K * std::exp(-M_PI * (a * a) * ((x * x) + (y * y)));
  float sinusoidal_carrier = std::cos(2.0 * M_PI * F_0 * ((x * std::cos(omega_0)) + (y * std::sin(omega_0))));
  return gaussian_envelop * sinusoidal_carrier;
}

unsigned morton(unsigned x, unsigned y)
{
  unsigned z = 0;
  for (unsigned i = 0; i < (sizeof(unsigned) * CHAR_BIT); ++i) {
    z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
  }
  return z;
}

class noise
{
public:
  noise(float K, float a, float F_0, float omega_0, float number_of_impulses_per_kernel, unsigned period, unsigned random_offset)
    :  K_(K), a_(a), F_0_(F_0), omega_0_(omega_0), period_(period), random_offset_(random_offset)
  {
    kernel_radius_ = std::sqrt(-std::log(0.05) / M_PI) / a_;
    impulse_density_ = number_of_impulses_per_kernel / (M_PI * kernel_radius_ * kernel_radius_);
  }
  float operator()(float x, float y) const
  {
    x /= kernel_radius_, y /= kernel_radius_;
    float int_x = std::floor(x), int_y = std::floor(y);
    float frac_x = x - int_x, frac_y = y - int_y;
    int i = int(int_x), j = int(int_y);
    float noise = 0.0;
    for (int di = -1; di <= +1; ++di) {
      for (int dj = -1; dj <= +1; ++dj) {
        noise += cell(i + di, j + dj, frac_x - di, frac_y - dj);
      }
    }
    return noise;
  }
  float cell(int i, int j, float x, float y) const
  {
    //unsigned s = (((unsigned(j) % period_) * period_) + (unsigned(i) % period_)) + random_offset_; // periodic noise
    unsigned s = morton(i, j); // nonperiodic noise
    if (s == 0) s = 1;
    pseudo_random_number_generator prng;
    prng.seed(s);
    double number_of_impulses_per_cell = impulse_density_ * kernel_radius_ * kernel_radius_;
    unsigned number_of_impulses = prng.poisson(number_of_impulses_per_cell);
    float noise = 0.0;

    
    //std::cout << i << ' ' << j << ' ' << x << ' ' << y << ' ' << s << ' ' << number_of_impulses << std::endl;
    for (unsigned i = 0; i < number_of_impulses; ++i) {
      float x_i = prng.uniform_0_1();
      float y_i = prng.uniform_0_1();
      float w_i = prng.uniform(-1.0, +1.0);
      float omega_0_i = prng.uniform(0.0, 2.0 * M_PI);
      float x_i_x = x - x_i;
      float y_i_y = y - y_i;
      //std::cout << "  " << x_i_x << ' ' << y_i_y << std::endl;
      if (((x_i_x * x_i_x) + (y_i_y * y_i_y)) < 1.0) {
        noise += w_i * gabor(K_, a_, F_0_, omega_0_, x_i_x * kernel_radius_, y_i_y * kernel_radius_); // anisotropic
        //std::cout << "     " << noise << std::endl;
        //noise += w_i * gabor(K_, a_, F_0_, omega_0_i, x_i_x * kernel_radius_, y_i_y * kernel_radius_); // isotropic
      }
    }
    return noise;
  }
  float variance() const
  {
    float integral_gabor_filter_squared = ((K_ * K_) / (4.0 * a_ * a_)) * (1.0 + std::exp(-(2.0 * M_PI * F_0_ * F_0_) / (a_ * a_)));
    return impulse_density_ * (1.0 / 3.0) * integral_gabor_filter_squared;
  }
private:
  float K_;
  float a_;
  float F_0_;
  float omega_0_;
  float kernel_radius_;
  float impulse_density_;
  unsigned period_;
  unsigned random_offset_;
};

// -----------------------------------------------------------------------------

#include <ctime>
#include <fstream>
#include <iostream>

int main(int argc, char* argv)
{
  std::cout << "Copyright (c) 2009 by Ares Lagae, Sylvain Lefebvre, George Drettakis and Philip Dutre." << std::endl;
  std::cout << "This code is designed for brevity and clarity. It is significantly slower than the GPU version." << std::endl;

  // ---------------------------------------------------------------------------

  unsigned resolution = 256;
  float K_ = 1.0;
  float a_ = 0.05;
  float F_0_ = 0.0625;
  float omega_0_ = M_PI / 4.0;
  float number_of_impulses_per_kernel = 64.0;
  unsigned period = 256;
  unsigned random_offset = std::time(0);

  // ---------------------------------------------------------------------------

  noise noise_(K_, a_, F_0_, omega_0_, number_of_impulses_per_kernel, period, random_offset);
  float* image = new float[resolution * resolution];
  float scale = 3.0 * std::sqrt(noise_.variance());
  for (unsigned row = 0; row < resolution; ++row) {
    for (unsigned col = 0; col < resolution; ++col) {
      float x = (float(col) + 0.5) - (float(resolution) / 2.0);
      float y = (float(row) + 0.5) - (float(resolution) / 2.0);

      float t = noise_(x, y);
      image[(row * resolution) + col] = 0.5 + (0.5 * (noise_(x, y) / scale));
    }
  }

  // ---------------------------------------------------------------------------

  std::ofstream file;
  file.open("noise.ppm", std::ios_base::binary);
  file << "P6" << '\n';
  file << resolution << ' ' << resolution << '\n';
  file << "255" << '\n';
  for (unsigned row = 0; row < resolution; ++row) {
    for (unsigned col = 0; col < resolution; ++col) {
      float intensity = image[(row * resolution) + col];
      unsigned char pixel[3];
      if (intensity <= 0.0) {
        pixel[0] = 0;
      }
      else if (intensity < 1.0) {
        pixel[0] = static_cast<unsigned char>(intensity * 255);
      }
      else {
        pixel[0] = 255;
      }
      pixel[1] = pixel[0];
      pixel[2] = pixel[0];
      file.write(reinterpret_cast<const char*>(&pixel), 3);
    }
  }

  // ---------------------------------------------------------------------------

  std::cout << "Please see 'noise.ppm'." << std::endl;
}
