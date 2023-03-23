//
// Created by kobayasi on 22-11-9.
//

#ifndef BAER_ETHERCAT_LEG_H
#define BAER_ETHERCAT_LEG_H

#include "calculate.h"

enum right_left {
    right = 0,
    left = 1
};

struct Length_Angle{
    double  OR;
    double  RL;
    double  LJ;
    double  JM;
    double  MK;
    double  KN;
    double  KB;
    double  BD;
    double  DG;
    double  GH;
    double  CDG;
    double  FGH;

    Length_Angle();
};

struct unit_link {
    int     ID;
    std::string     name;
    std::string     point;
    Vec3<double>    p;
    Mat3<double>    R;
    Mat3<double>    R0;
    Vec3<double>    v;
    Vec3<double>    w;
    double          q;
    double          dq;
    double          ddq;
    Vec3<double>    a;
    Vec3<double>    b;
    Vec3<double>    vertex;
    double          face;
    double          m;
    Vec3<double>    c;
    Vec3<double>    I;
};

struct TerminalClass {
    double  x;
    double  y;
    double  z;
    double  p;
    double  yaw_ik;
    double  roll_ik;

    double theta[5];
};

class LegClass {
public:
    LegClass();
    void FKinematics(double* theta);
    void IKinematics(double* Tra, double* pas);

    double Terminal[2][5];
    double theta[2][7];
protected:
    Length_Angle    A;
    unit_link   link[2][7];
};




#endif //BAER_ETHERCAT_LEG_H
