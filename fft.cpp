#include "gabor.hpp"

const int maxn = 1024;
Complex w[maxn];
void fft(vector<Complex> &y, int bit, int on = 1) {
  for (int i = 0, j = 0; i < (1 << bit); i++) {
    if (i > j) swap(y[i], y[j]);
    for (int l = (1 << (bit - 1)); (j ^= l) < l; l >>= 1);
  }

	int n = (1 << bit);
  w[0] = Complex(1, 0);
	for (int s = 0; s < bit; s++) {
		int m = (1 << s) << 1, p = (1 << s);
		Complex g = Complex(cos(on * 2 * M_PI / m), sin(on * 2 * M_PI / m));
		for(int j = p; j >= 0; j -= 2) w[j] = w[j >> 1];
		for(int j = 1; j < p; j += 2) w[j] = w[j - 1] * g;

		for (int i = 0; i < n; i += m) {
      Complex *a = &y[i], *b = &y[i + p];
			for (int j = 0; j < p; j++, a++, b++) {
				Complex s = *b * w[j];
				*b = *a - s;
				*a = *a + s;
			}
		}
	}

	if (on == -1) 
		for (int i = 0; i < n; i++) y[i] /= n;
}

class Poly {
public:
	int n;

	vector<Complex> a;
	Poly() {}
	Poly(int _n, vector<Complex> data) {
		// less than
		// n must be a power of 2
		n = _n; 

		int bit = 0;
		while ((1 << bit) < n) bit++;

		a.swap(data);
		fft(a, bit);
	}

	Complex operator[] (const int &idx) {
		return a[idx];
	}
} ;

int _n;
float_type adjust(int x) {
	return (float_type)x - _n / 2;
}

vector< vector<float_type> > twodimsfft(int n, vector< vector<float_type> > data) {
	_n = n;
	vector< vector<float_type> > rt(n);
	for (int i = 0; i < n; i++) rt[i].resize(n);

	for (int i = 0; i < n / 2; i++) 
	for (int j = 0; j < n; j++) swap(data[i][j], data[i + n / 2][j]);
	for (int i = 0; i < n; i++) 
	for (int j = 0; j < n / 2; j++) swap(data[i][j], data[i][j + n / 2]);


	// printf("Start two dimensions fft!\n");
	vector<Poly> polys(n);
	vector<Complex> _t(n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) _t[j] = Complex(data[i][j], 0);
		polys[i] = Poly(n, _t);
	}
	// printf("Finished Calculating polys!\n");

	vector<Complex> tmp(n);
	Poly _poly;
	for (int v = 0; v < n; v++) {
		for (int x = 0; x < n; x++) tmp[x] = polys[x][v];
		_poly = Poly(n, tmp);
		for (int u = 0; u < n; u++) 
			rt[u][v] = _poly[u].real();
	}
	// printf("Finished FFT!\n");

	for (int i = 0; i < n / 2; i++) 
	for (int j = 0; j < n; j++) swap(rt[i][j], rt[i + n / 2][j]);
	for (int i = 0; i < n; i++) 
	for (int j = 0; j < n / 2; j++) swap(rt[i][j], rt[i][j + n / 2]);

	return rt;
}

