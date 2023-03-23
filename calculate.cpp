//
// Created by kobayasi on 22-10-29.
#include "calculate.h"
#include "cmath"

double Sign(double x){
    if (x < 0)
        return -1.0;
    else if (x > 0)
        return 1.0;
    else
        return 0.0;
}

double Ramp(double dst, double v0, double inc){
    if (fabs(v0 - dst) < inc)
        return v0;
    else
        return (v0 + Sign(dst - v0)*inc);
}

double deg2rad(double x)
{
    double a;
    a = x/180.0*PI;
    return a;
}

double rad2deg(double x)
{
    double a;
    a = x*180.0/PI;
    return a;
}

double limit_cal(double x, double max, double min){
    if (x < min)        return min;
    else if (x > max)    return max;
    else                    return x;
}

void v1p2::cal(int index) {
    double con;
    con = b*b - 4*a*c;
    if (con > 0.0000000001 | con < -0.0000000001){
        ero = false;
        if (index == 0)
            sol = (-b + sqrt(con))/(2*a);
        else if (index == 1)
            sol = (-b - sqrt(con))/(2*a);
        else
            sol = 0;
    }
    else if (con < 0.0000000001 & con > 0.0000000001){
        ero = false;
        sol = -b/(2*a);
    }
    else{
        ero = true;
        sol = 0;
    }
}

void v1p2_linear::cal(int index) {
    double a, b, c;
    double x[2], y[2];
    double con;
    a = (1+k*k);
    b = 2*k*t;
    c = t*t-h*h;
    con = b*b - 4*a*c;
    if (con < 0)
        ero = true;
    else
        ero = false;
    x[0] = (-b + sqrt(con))/(2*a);
    x[1] = (-b - sqrt(con))/(2*a);

    y[0] = k*x[0] + t;
    y[1] = k*x[1] + t;

    output[0] = x[index];
    output[1] = y[index];
}

Mat3<double> Rx3(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Mat3<double> M;
    M << 1, 0, 0,
            0, cos(theta), -sin(theta),
            0, sin(theta),  cos(theta);
    return M;
}

Mat3<double> Ry3(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Mat3<double> M;
    M << cos(theta), 0, sin(theta),
            0,             1,               0,
            -sin(theta), 0, cos(theta);
    return M;
}

Mat3<double> Rz3(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Mat3<double> M;
    M << cos(theta), -sin(theta), 0,
            sin(theta), cos(theta),  0,
            0,             0,              1;
    return M;
}

Mat2<double> Rz(double theta){
    // for 2D-XY vector, rotation matrix along z axis
    Mat2<double> M;
    M << cos(theta),-sin(theta),
            sin(theta), cos(theta);
    return M;

}

Mat3<double> cross_product(Vec3<double> a){
    Mat3<double> aa;
    aa << 0, -a(2), a(1),
            a(2), 0, -a(0),
            -a(1), a(0), 0;

    return aa;
}

Mat3<double> Rodrigues(Vec3<double> a, double theta){
    Mat3<double> R;
    Mat3<double> E;

    E << 1, 0, 0, 0, 1, 0, 0, 0, 1;
    R = E + cross_product(a)*sin(theta) + cross_product(a)*cross_product(a)*(1-cos(theta));

    return R;
}

double approach(double ini, double dst, double s){
    double aa;
    if (s < 1)
        aa = ini + (dst - ini)*(0.5 - 0.5*cos(2*PI*s/2));
    else
        aa = dst;

    return aa;
}

void filter::LPF_cal(double fre) {
    double  Og = 2;
    double  w = 2*PI*fre;
    double  LPF_para1, LPF_para2, LPF_para3;
    LPF_para1 = w*w*Ts*Ts/(1 + Og*w*Ts + w*w*Ts*Ts);
    LPF_para2 = (2 + Og*w*Ts)/(1 + Og*w*Ts + w*w*Ts*Ts);
    LPF_para3 = -1/(1 + Og*w*Ts + w*w*Ts*Ts);

    out = LPF_para1 * in + LPF_para2 * out_pre[0] + LPF_para3 * out_pre[1];

    in_pre[2] = in_pre[1];
    in_pre[1] = in_pre[0];
    in_pre[0] = in;
    out_pre[2] = out_pre[1];
    out_pre[1] = out_pre[0];
    out_pre[0] = out;

}

void filter::notch_cal(double fre){
    double  num[3], den[3];
    double  a[3], b[3];
    double  w = 2*PI*fre;
    double  k1, k2;

    a[0] = w*w*Ts*Ts + 4*Ts*k1*w + 4;
    a[1] = 2*w*w*Ts*Ts - 8;
    a[2] = w*w*Ts*Ts - 4*Ts*k1*w + 4;

    b[0] = w*w*Ts*Ts + 4*Ts*k2*w + 4;
    b[1] = 2*w*w*Ts*Ts - 8;
    b[2] = w*w*Ts*Ts - 4*Ts*k2*w + 4;

    out = num[0]*in + num[1]*in_pre[0] + num[2]*in_pre[1] -(den[1]*out_pre[0] + den[2]*out_pre[0]);

    in_pre[2] = in_pre[1];
    in_pre[1] = in_pre[0];
    in_pre[0] = in;
    out_pre[2] = out_pre[1];
    out_pre[1] = out_pre[0];
    out_pre[0] = out;
}

void PID_STRUCT::cal() {
    if (enable == 1){
        err = ref - fpb;

        //P calculate
        pval = sign*Kp*err;

        //I calculate
        Ierr = err;
        if (Ierr > dead_zone[0])
            Ierr = Ierr - dead_zone[0];
        else if (Ierr < dead_zone[1])
            Ierr = Ierr - dead_zone[1];
        else
            Ierr = 0.0;
        ival = ival + sign*Ki*Ierr*Ts;
        ival = limit_cal(ival, limit[0], limit[1]);

        //D calculate
        if (err_pre == 0) // prevent the impact at initial time
            err_pre = err;
        //if (err == 0) // prevent dval false caused by the same data
        //    err = err_pre;
        dval_LPF.in = sign*Kd*(err - err_pre)/Ts;
        dval_LPF.LPF_cal(dval_fre);
        dval = dval_LPF.out;

        //val calculate
        val_notch.in = pval + ival + dval;
        val_notch.notch_cal(val_fre);
        val = limit_cal(val_notch.out, limit[0], limit[1]);


        err_pre = err;
        dval_pre = dval;
    }
    else{
        val = 0.0;
        pval = 0.0;
        ival = 0.0;
        dval = 0.0;
        err_pre = 0.0;
        dval_LPF.in = 0.0;
        dval_LPF.out = 0.0;
        for (int i = 0; i < 3; i++){
            dval_LPF.in_pre[i] = 0.0;
            dval_LPF.out_pre[i] = 0.0;
        }
    }
}

//
