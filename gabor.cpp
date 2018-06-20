typedef double float_type;
class pseudo_random_number_generator {
private:
    unsigned x;
public:
    void seed(unsigned s) {x = s;}
    unsigned next() {return x *= 3039177861u;}
    float_type uniform_0_1() {return float_type(next()) / float_type(UINT_MAX);}
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
} ;

class GaborNoise {
public:
    static constexpr int bit = 1;
    static constexpr int blocks = (1 << bit);
    //static constexpr float_type kernel_radius = 50;
    //static constexpr float_type kernel_radius = float_type(n) / blocks;
    static constexpr unsigned impluse_density_ = 64;
    static constexpr float_type K = 1.0;
    static constexpr float_type a = 0.05;
    static constexpr float_type F_0 = .0625;
    static constexpr float_type lower_limit = 0;
    static constexpr float_type upper_limit = 2 * M_PI;

    static constexpr float_type kernel_radius = sqrt(-log(0.05) / M_PI) / a;
    static constexpr unsigned number_of_impulses_per_cell = 64.0;

    static constexpr float_type omega_0 = M_PI / 4.0;

public:
    float_type frac(float_type x) {return x - floor(x);}
    inline unsigned morton(unsigned x, unsigned y) {
        unsigned z = 0;
        for (unsigned i = 0; i < (sizeof(unsigned) * CHAR_BIT); ++i) {
            z |= ((x & (1 << i)) << i) | ((y & (1 << i)) << (i + 1));
        }
        return z;
    }

    inline float_type sqr(float_type x) {return x * x;}

    float_type gaussian(float_type x, float_type y, float_type K, float_type a) {
        return K * exp(-M_PI * sqr(a) * (sqr(x) + sqr(y)));
    }

    inline float_type gabor(float_type K, float_type a, float_type F_0, float_type omega_0, float_type x, float_type y) {
    	return K * exp(-M_PI * sqr(a) * (sqr(x) + sqr(y))) * cos(2 * M_PI * F_0 * (x * cos(omega_0) + y * sin(omega_0)));
    }

    typedef pair<unsigned, unsigned> Point;
    map<Point, float_type> mp;
    float_type cell(int i, int j, float_type x, float_type y) {

        unsigned s = morton(i, j);
        if (s == 0) s = 1;
        pseudo_random_number_generator prng;
        prng.seed(s);
        unsigned number_of_impulses = prng.poisson(number_of_impulses_per_cell / M_PI);

        float_type sum = 0.0;
        float_type x_i, y_i, w_i, omega_0_i;
        for (unsigned i = 0; i < number_of_impulses; i++) { 
            x_i = prng.uniform_0_1(), y_i = prng.uniform_0_1();
            w_i = prng.uniform(-1.0, +1.0);
            omega_0_i = prng.uniform(lower_limit, upper_limit);

            float_type x_i_x = x - x_i;
            float_type y_i_y = y - y_i;
            //cout << "  " << x_i_x << ' ' << y_i_y << endl;
            if (((x_i_x * x_i_x) + (y_i_y * y_i_y)) < 1.0) {
                sum += w_i * gabor(K, a, F_0, omega_0, (x - x_i) * kernel_radius, (y - y_i) * kernel_radius);
                //cout << "     " << sum << endl;
                //sum += w_i * gabor(K, a, F_0, omega_0_i, (x - x_i) * kernel_radius, (y - y_i) * kernel_radius); 
            }
        }
        return sum;
    }
public:
    float_type noise(float_type x, float_type y) {
        x /= kernel_radius, y /= kernel_radius;
        float_type sum = 0.0;
        for (int i = -1; i <= 1; i++) 
            for (int j = -1; j <= 1; j++) 
                sum += cell(int(floor(x)) + i, int(floor(y)) + j, frac(x) - i, frac(y) - j);
        least = min(least, sum);
        best = max(best, sum);
        return sum;
    }

    float_type least, best;
    unsigned random_offset;
} ;