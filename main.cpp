#include <bits/stdc++.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include "gabor.hpp"
using namespace cv;
using namespace std;

const int n = 512;

Mat src(2 * n, 2 * n, CV_8UC1);
// GaborNoise_isotropic gn;
GaborNoise_anisotropic gn_ani;
GaborNoise_isotropic gn_iso;
GaborNoise_isotropic_kernel gn_iso_kernel;

GaborNoise *gn = &gn_iso;

typedef vector<vector<float_type> > LIST2D;

LIST2D img[4];

void calc(int stx, int edx, int sty, int edy) {
    for (int i = stx; i < edx; i++) 
    for (int j = sty; j < edy; j++) 
        img[0][i][j] = gn->noise(.5 + i, .5 + j);
}

void draw(const LIST2D &tmp, int x, int y) {
    float_type least = 10000000, best = -10000000;
    for (int i = 0; i < n; i++) 
    for (int j = 0; j < n; j++) { 
        least = min(least, tmp[i][j]);
        best = max(best, tmp[i][j]);
    }

    for (int i = 0; i < n; i++) {
        uchar *p = src.ptr<uchar>(i + x);
        for (int j = 0; j < y; j++, p++);
        for (int j = 0; j < n; j++, p++) 
            *p = static_cast<uchar>((tmp[i][j] - least) / (best - least) * 255.0);
    }
}

void draw() {
    int st = clock();
    thread t[4];
    t[0] = thread(calc, 0, n / 2, 0, n / 2);
    t[1] = thread(calc, 0, n / 2, n / 2, n);
    t[2] = thread(calc, n / 2, n, 0, n / 2);
    t[3] = thread(calc, n / 2, n, n / 2, n);
    for (int i = 0; i < 4; i++) t[i].join();

    for (int i = 0; i < n; i++) 
    for (int j = 0; j < n; j++) 
        img[1][i][j] = gn->get_gabor(5 * (i - n / 2.0) / n, 5 * (j - n / 2.0) / n);
    
    draw(img[0], 0, 0);
    draw(img[1], 0, n);
    img[2] = twodimsfft(n, img[0]);
    draw(img[2], n, 0);



    for (int i = 0; i < n; i++) 
    for (int j = 0; j < n; j++) 
        img[3][i][j] = gn->fourier((i - n / 2.0) / n, (j - n / 2.0) / n); 
    
    draw(img[3], n, n);
    cout << (clock() - st) / 1000.0 << endl;
}

int flag;
int _K = 100, _a = 5, _F_0 = 625, _number_of_impulses_per_cell = 64, _error_ratio = 5;
int _omega = 0;
int _lower_bound, _upper_bound = 100;

void onchange(int _, void *__) {
    gn->set(_K / 100.0, _a / 100.0, _F_0 / 10000.0, _number_of_impulses_per_cell, _error_ratio / 100.0);
    gn->set_bound(_lower_bound / 100.0 * (2 * M_PI), _upper_bound / 100.0 * (2 * M_PI));
    // if (gn->set_omega) {
    //     gn->set_omega((_omega / 100.0) * (2 * M_PI));
    // }
    //draw();
    //gn.set_F_0(_F_0 / 10000.0);
}

void omega_changed(int _, void *__) {
    gn_ani.set_omega((_omega / 100.0) * (M_PI));
}

void submit(int _, void *__) {
    if (_ == 0) {
        gn = &gn_ani;
        printf("using anisotropic.\n");
    }
    else if (_ == 1) {
        gn = &gn_iso;
        printf("using isotropic.\n");
    }
    else if (_ == 2) {
        gn = &gn_iso_kernel;
        printf("using isotropic with isotropic kernel.\n");
    }
}

int main() {
    namedWindow("demo");

    createTrackbar("K", "demo", &_K, 100, onchange, NULL);
    createTrackbar("a", "demo", &_a, 20, onchange, NULL);
    createTrackbar("F_0", "demo", &_F_0, 1000, onchange, NULL);
    createTrackbar("npc", "demo", &_number_of_impulses_per_cell, 100, onchange, NULL);
    createTrackbar("error", "demo", &_error_ratio, 10, onchange, NULL);
    createTrackbar("lower", "demo", &_lower_bound, 100, onchange, NULL);
    createTrackbar("upper", "demo", &_upper_bound, 100, onchange, NULL);
    createTrackbar("omega", "demo", &_omega, 100, omega_changed, NULL);
    
    flag = 1;
    createTrackbar("submit", "demo", &flag, 2, submit, NULL);
    
    for (int i = 0; i < 4; i++) {
        img[i].resize(n);
        for (int j = 0; j < n; j++) img[i][j].resize(n);
    }

    while (true) {
        draw();
        imshow("demo", src);
        char ch = waitKey();
        if (ch == 's') {
            imwrite("noise.png", src(Rect(0, 0, n, n)));
        }
    }
    return 0;
}
