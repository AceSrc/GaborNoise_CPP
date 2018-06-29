#include <bits/stdc++.h>
using namespace std;

typedef double float_type;
typedef complex<float_type> Complex;

vector<vector<float_type> > twodimsfft(int n, vector<vector<float_type> > data);

typedef double float_type;
class pseudo_random_number_generator {
 private:
  unsigned x;

 public:
  void seed(unsigned s) { x = s; }
  unsigned next() { return x *= 3039177861u; }
  float_type uniform_0_1() { return float_type(next()) / float_type(UINT_MAX); }
  float_type uniform(float_type min, float_type max) {
    return min + uniform_0_1() * (max - min);
  }
  unsigned poisson(float_type mean) {
    float_type g = exp(-mean);
    unsigned rt = 0;
    float_type t = uniform_0_1();
    while (t > g) {
      rt++;
      t *= uniform_0_1();
    }
    return rt;
  }
};

class GaborNoise {
 protected:
  // unsigned impluse_density;
  float_type K;
  float_type a;
  float_type F_0;
  unsigned number_of_impulses_per_cell;
  float_type error_ratio;

  float_type lower_bound, upper_bound;

  float_type kernel_radius;
  unsigned random_offset;

  inline float_type frac(float_type x) { return x - floor(x); }
  inline unsigned morton(unsigned x, unsigned y) {
    unsigned z = 0;
    for (unsigned i = 0; i < (sizeof(unsigned) * CHAR_BIT); ++i) {
      z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
    }
    return z;
  }

  inline float_type sqr(float_type x) const { return x * x; }

  /*float_type gaussian(float_type x, float_type y, float_type K, float_type a)
  { return K * exp(-M_PI * sqr(a) * (sqr(x) + sqr(y)));
  }*/

  virtual float_type gabor(float_type K, float_type a, float_type F_0,
                           float_type omega_0, float_type x,
                           float_type y) const {
    return 0;
  }

  virtual float_type variance() const { return 0; }

  /*{
      return K * exp(-M_PI * sqr(a) * (sqr(x) + sqr(y))) * cos(2 * M_PI * F_0 *
  (x * cos(omega_0) + y * sin(omega_0)));
  }*/

  virtual float_type fourier(float_type x, float_type y) const { return 0; }
 public:
  void set(float_type _K = 1.0, float_type _a = .05, float_type _F_0 = .0625,
           unsigned _number_of_impulses_per_cell = 64,
           float_type _error_ratio = .05) {
    K = _K;
    a = _a;
    F_0 = _F_0;
    number_of_impulses_per_cell = _number_of_impulses_per_cell;
    error_ratio = _error_ratio;
  }

  void set_a(float_type _a = .05) { a = _a; }

  void set_F_0(float_type _F_0 = .0625) { F_0 = _F_0; }

  void set_bound(float_type _lower_bound, float_type _upper_bound) {
    lower_bound = _lower_bound;
    upper_bound = _upper_bound;
  }

  GaborNoise(float_type _K = 1.0, float_type _a = .05, float_type _F_0 = .0625,
             unsigned _number_of_impulses_per_cell = 64,
             float_type _error_ratio = .05)
      : K(_K),
        a(_a),
        F_0(_F_0),
        number_of_impulses_per_cell(_number_of_impulses_per_cell),
        error_ratio(_error_ratio) {
    random_offset = rand();
    lower_bound = 0.0;
    upper_bound = 2 * M_PI;
  }

  float_type cell(int i, int j, float_type x, float_type y) {
    unsigned s = morton(i, j) + random_offset;
    if (s == 0) s = 1;
    pseudo_random_number_generator prng;
    prng.seed(s);
    unsigned number_of_impulses =
        prng.poisson(number_of_impulses_per_cell / M_PI);

    float_type sum = 0.0;
    float_type x_i, y_i, w_i, omega_0_i;
    for (unsigned i = 0; i < number_of_impulses; i++) {
      x_i = prng.uniform_0_1(), y_i = prng.uniform_0_1();
      w_i = prng.uniform(-1.0, +1.0);
      omega_0_i = prng.uniform(lower_bound, upper_bound);

      float_type x_i_x = x - x_i;
      float_type y_i_y = y - y_i;
      if (((x_i_x * x_i_x) + (y_i_y * y_i_y)) < 1.0) {
        sum += w_i * gabor(K, a, F_0, omega_0_i, (x - x_i) * kernel_radius,
                           (y - y_i) * kernel_radius);
      }
    }
    return sum;
  }

  float_type get_kernel_radius() const {
    return sqrt(-log(error_ratio) / M_PI) / a;
  }

  float_type get_kernel(float_type x, float_type y) const {
    return gabor(K, a, F_0, lower_bound, x * kernel_radius, y * kernel_radius);
  }

  float_type get_fourier(float_type x, float_type y) const {
    return fourier(x, y);
  }

  void set_error_ratio(float_type _error_ratio) {
    error_ratio = _error_ratio;
    kernel_radius = sqrt(-log(error_ratio) / M_PI) / a;
  }

  float_type noise(float_type x, float_type y) {
    kernel_radius = sqrt(-log(error_ratio) / M_PI) / a;
    x /= kernel_radius, y /= kernel_radius;

    // cout << "H" << endl;
    // cout << K << ' ' << a << ' ' << F_0 << ' ' << number_of_impulses_per_cell
    // << ' ' << kernel_radius << endl;
    float_type sum = 0.0;
    for (int i = -1; i <= 1; i++)
      for (int j = -1; j <= 1; j++)
        sum += cell(int(floor(x)) + i, int(floor(y)) + j, frac(x) - i,
                    frac(y) - j);
    return sum;
  }
};

class GaborNoise_isotropic : public GaborNoise {
 private:
  float_type gabor(float_type K, float_type a, float_type F_0,
                   float_type omega_0, float_type x, float_type y) const {
    return K * exp(-M_PI * sqr(a) * (sqr(x) + sqr(y))) *
           cos(2 * M_PI * F_0 * (x * cos(omega_0) + y * sin(omega_0)));
  }

 public:
  virtual float_type fourier(float_type x, float_type y) const {
    float_type f_r = sqrt(sqr(x) + sqr(y));
    return sqr(K) / (4 * sqrt(2) * M_PI * F_0 * sqr(a) * a) *
           exp(-2 * M_PI / sqr(a) * sqr(f_r - F_0));
  }

  float_type variance() const { return 0; }
};

class GaborNoise_anisotropic : public GaborNoise {
 private:
  float_type omega;
  float_type gabor(float_type K, float_type a, float_type F_0,
                   float_type omega_0, float_type x, float_type y) const {
    return K * exp(-M_PI * sqr(a) * (sqr(x) + sqr(y))) *
           cos(2 * M_PI * F_0 * (x * cos(omega) + y * sin(omega)));
  }

  virtual float_type fourier(float_type x, float_type y) const {
    return K / (2 * sqr(a)) *
           (exp(-M_PI / sqr(a) *
                (sqr(x - F_0 * cos(omega)) + sqr(y - F_0 * sin(omega)))) +
            exp(-M_PI / sqr(a) *
                (sqr(x + F_0 * cos(omega)) + sqr(y + F_0 * sin(omega)))));
  }

 public:
  void set_omega(float_type _omega = M_PI / 4.0) { this->omega = _omega; }

  float_type variance() const { return 0; }

  float_type get_kernel(float_type x, float_type y) const {
    return gabor(K, a, F_0, omega, x * kernel_radius, y * kernel_radius);
  }
};

class GaborNoise_isotropic_kernel : public GaborNoise {
 private:
  float_type A, B, C;
  float_type I(float_type x) const {
    float_type ax, ans;
    float_type y;
    if ((ax = fabs(x)) < 3.75) {
      y = x / 3.75;
      y *= y;
      ans =
          1.0 +
          y * (3.5156229 +
               y * (3.0899424 +
                    y * (1.2067492 + y * (0.2659732 + y * (0.360768e-1 +
                                                           y * 0.45813e-2)))));
    } else {
      y = 3.75 / ax;
      ans = (exp(ax) / sqrt(ax)) *
            (0.39894228 +
             y * (0.1328592e-1 +
                  y * (0.225319e-2 +
                       y * (-0.157565e-2 +
                            y * (0.916281e-2 +
                                 y * (-0.2057706e-1 +
                                      y * (0.2635537e-1 +
                                           y * (-0.1647633e-1 +
                                                y * 0.392377e-2))))))));
    }
    return ans;
  }

  float_type J(float_type x) const {
    float ax, z;
    double xx, y, ans, ans1, ans2;
    if ((ax = fabs(x)) < 8.0) {
      y = x * x;
      ans1 = 57568490574.0 +
             y * (-13362590354.0 +
                  y * (651619640.7 +
                       y * (-11214424.18 +
                            y * (77392.33017 + y * (-184.9052456)))));
      ans2 = 57568490411.0 +
             y * (1029532985.0 +
                  y * (9494680.718 +
                       y * (59272.64853 + y * (267.8532712 + y * 1.0))));
      ans = ans1 / ans2;
    } else {
      z = 8.0 / ax;
      y = z * z;
      xx = ax - 0.785398164;
      ans1 = 1.0 + y * (-0.1098628627e-2 +
                        y * (0.2734510407e-4 +
                             y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
      ans2 = -0.1562499995e-1 +
             y * (0.1430488765e-3 +
                  y * (-0.6911147651e-5 +
                       y * (0.7621095161e-6 - y * 0.934945152e-7)));
      ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
    }
    return ans;
  }

  float_type gabor(float_type K, float_type a, float_type F_0,
                   float_type omega_0, float_type x, float_type y) const {
    float_type r = sqrt(sqr(x) + sqr(y));

    // cout << K << ' ' << exp(A * sqr(r)) << ' ' << B << ' ' << J(C * r) <<
    // endl;

    // return K * exp(A * sqr(r)) / (sqr(a));
    // return 2 * M_PI * J(2 * M_PI * F_0 * r) * F_0;
    return K * exp(A * sqr(r)) * B * J(C * r);
  }


  virtual float_type fourier(float_type x, float_type y) const {
    float_type f_r = sqrt(sqr(x) + sqr(y));
    float_type A = 2 * M_PI * K * F_0 / sqr(a);
    float_type B = -M_PI / sqr(a) * (sqr(f_r) + sqr(F_0));
    float_type C = 2 * M_PI * F_0 / sqr(a) * f_r;
    return A * exp(B) * I(C);
  }

 public:
  GaborNoise_isotropic_kernel(float_type _K = 1.0, float_type _a = .05,
                              float_type _F_0 = .0625,
                              unsigned _number_of_impulses_per_cell = 64,
                              float_type _error_ratio = .05)
      : GaborNoise(_K, _a, _F_0, _number_of_impulses_per_cell, _error_ratio) {
    A = -M_PI * sqr(a);
    B = 2 * M_PI * F_0;
    C = 2 * M_PI * F_0;
  }

  float_type variance() const { return 0; }
};
