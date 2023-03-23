//
// Created by kobayasi on 22-10-29.
//

#ifndef BAER_ETHERCAT_CALCULATE_H
#define BAER_ETHERCAT_CALCULATE_H

#define PI 3.14159265359

#include "cppTypes.h"

struct v1p2 {
    // a*x^2 + b*x + c = 0
    double  a;
    double  b;
    double  c;
    bool    ero;
    double  sol;

    void cal(int index);
};

struct v1p2_linear {
    // 2 : x^2 + y^2 = h^2
    // 1 : y = k*x + t
    double  h;
    double  k;
    double  t;
    bool    ero;
    double  output[2];

    void cal(int index);
};

struct filter {
    double  in;
    double  in_pre[3];
    double  out;
    double  out_pre[3];

    void LPF_cal(double fre);
    void notch_cal(double fre);
};

struct PID_STRUCT {
    double  Kp;
    double  Ki;
    double  Kd;
    double  limit[2];//up,down
    double  dead_zone[2];//up,down
    double  sign;
    double  dval_fre;
    double  val_fre;

    double  ref;
    double  fpb;
    double  err;
    double  Ierr;
    double  err_pre;
    double  dval_pre;
    double  val;
    double  pval;
    double  ival;
    double  dval;
    filter  dval_LPF;
    filter  val_notch;

    uint16_t enable;

    void cal();
};


extern double  Ts;

double Sign(double x);
double Ramp(double dst, double v0, double inc);
double deg2rad(double x);
double rad2deg(double x);
double limit_cal(double x, double max, double min);
Mat3<double> Rx3(double theta);
Mat3<double> Ry3(double theta);
Mat3<double> Rz3(double theta);
Mat2<double> Rz(double theta);
Mat3<double> Rodrigues(Vec3<double> a, double theta);
Mat3<double> cross_product(Vec3<double> a);

double approach(double ini, double dst, double s);

#endif //BAER_ETHERCAT_CALCULATE_H
