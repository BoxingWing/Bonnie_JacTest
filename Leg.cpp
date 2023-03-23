//
// Created by kobayasi on 22-11-9.
//
#include <iostream>
#include <fstream>

#include "Leg.h"

//14

Length_Angle::Length_Angle() {
    OR = 0.125;
    RL = 0.02;
    LJ = 0.0585;
    JM = 0.071;//0.065;//
    MK = 0.022;
    KN = 0.0557;//0.0565;//
    KB = 0.0647;//0.065;//
    BD = 0.168;//0.158;//
    DG = 0.23;
    GH = 0.2;

    CDG = 3.2812;//3.2563;
    FGH = 2.8623;//3.1915;
}

LegClass::LegClass(){
    std::string name[2][7] = {"roll_r", "yaw_r", "pitch_r", "knee_r", "pas1_r", "pas2_r", "ankle_r",
                              "roll_l", "yaw_l", "pitch_l", "knee_l", "pas1_l", "pas2_l", "ankle_l"};

    std::string point[2][7] = {"OJ", "JK", "KBD", "AB", "CDG", "FGH", "Hh",
                               "OJ", "JK", "KBD", "AB", "CDG", "FGH", "Hh"};
    for (int l = 0; l < 2; l++) {
        for (int i; i < 7; i++) {
            link[l][i].ID = i;
            link[l][i].point = point[l][i];
            link[l][i].name = name[l][i];
        }
        link[l][0].p << 0, 0, 0;

        link[l][0].R0 << 0, 0, -1, 0, -1, 0, -1, 0, 0;
        link[l][1].R0 << 0, 0, -1, 0, -1, 0, -1, 0, 0;
        if (l == 0)     link[l][2].R0 << 0, 1, 0, 0, 0, -1, -1, 0, 0;
        else            link[l][2].R0 << 0, -1, 0, 0, 0, 1, -1, 0, 0;
        if (l == 0)     link[l][3].R0 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
        else            link[l][3].R0 << 0, 1, 0, -1, 0, 0, 0, 0, 1;
        link[l][4].R0 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        link[l][5].R0 << 1, 0, 0, 0, 1, 0, 0, 0, 1;
        if (l == 0)     link[l][6].R0 << 0, -1, 0, 1, 0, 0, 0, 0, 1;
        else            link[l][6].R0 << 0, 1, 0, -1, 0, 0, 0, 0, 1;

        link[l][0].a << -1, 0, 0;
        link[l][1].a << -1, 0, 0;
        if (l == 0)     link[l][2].a << 0, -1, 0;
        else            link[l][2].a << 0, 1, 0;
        link[l][3].a << 0, 0, 1;
        link[l][4].a << 0, 0, 1;
        link[l][5].a << 0, 0, 1;
        link[l][6].a << 0, 0, 1;

//        if (l == 0)     link[l][0].b << 0, -A.OR, 0;
//        else    link[l][0].b << 0, A.OR, 0;
        link[l][0].b << 0, 0, 0;
        link[l][1].b << A.RL, 0, -A.LJ;
        if (l == 0)     link[l][2].b << 0, -A.MK, -A.JM;
        else            link[l][2].b << 0, A.MK, -A.JM;
        link[l][3].b << 0, 0, -A.KN;
        link[l][4].b << A.BD, 0, -A.KB;
        link[l][5].b << A.DG, 0, 0;
        link[l][6].b << A.GH, 0, 0;
    }
}

void LegClass::FKinematics(double* theta){
    int sign;
    Vec3<double> vfo;
    double tmp;
    link[0][0].q = theta[0];//deg2rad(12);
    link[0][1].q = theta[1];//deg2rad(20);
    link[0][2].q = theta[2];//deg2rad(-10);
    link[0][3].q = theta[3];//deg2rad(120);
    link[0][4].q = theta[5];//deg2rad(-3);
    link[0][5].q = theta[6];//deg2rad(5);
    link[0][6].q = theta[4];//deg2rad(90);
    link[1][0].q = theta[7];//deg2rad(12);
    link[1][1].q = theta[8];//deg2rad(20);
    link[1][2].q = theta[9];//deg2rad(-10);
    link[1][3].q = theta[10];//deg2rad(120);
    link[1][4].q = theta[12];//deg2rad(-3);
    link[1][5].q = theta[13];//deg2rad(5);
    link[1][6].q = theta[11];//deg2rad(90);
    for(int l = 0; l < 2; l ++){
        sign = (int)(1 - l*2);

        link[l][0].p = link[l][0].b;
//        printf("p0 = %.4f, %.4f, %.4f \n", link[l][0].p(0), link[l][0].p(1), link[l][0].p(2));
        link[l][0].R = link[l][0].R0*Rz3(link[l][0].q);

        link[l][1].p = link[l][0].p + link[l][0].R*link[l][1].b;//(A.OJ, 0, 0);
//        printf("p1 = %.4f, %.4f, %.4f \n", link[l][1].p(0), link[l][1].p(1), link[l][1].p(2));
        link[l][1].R = link[l][0].R*link[l][1].R0*Rz3(link[l][1].q);//Vec3<double> Vec_roll = Rx3(-1*roll)*Vec3<double>(0, 0, -1);

        link[l][2].p = link[l][1].p + link[l][1].R*link[l][2].b;//Vec3<double> PK = PJ + A.JK*Vec_roll;
//        printf("p2 = %.4f, %.4f, %.4f \n", link[l][2].p(0), link[l][2].p(1), link[l][2].p(2));
        link[l][2].R = link[l][1].R*link[l][2].R0*Rz3(link[l][2].q);//Vec3<double> Vec_yaw = Rx3(-1*roll)*Rz3(-1*yaw)*Vec3<double>(0, 1, 0);

        link[l][3].p = link[l][2].p + link[l][2].R*link[l][3].b;//Vec3<double> PB = PK + A.KB*Vec_yaw;
//        printf("p3 = %.4f, %.4f, %.4f \n", link[l][3].p(0), link[l][3].p(1), link[l][3].p(2));
        link[l][3].R = link[l][2].R*link[l][3].R0*Rz3(link[l][3].q);//Vec3<double> va1 = Ry3(-1*theta1)*Vec3<double>(0, 0, -1);Vec3<double> vb1 = Ry3(-1*theta2)*va1;

        link[l][4].p = link[l][2].p + link[l][2].R*link[l][4].b;
//        printf("p4 = %.4f, %.4f, %.4f \n", link[l][4].p(0), link[l][4].p(1), link[l][4].p(2));
        link[l][4].R = link[l][2].R*link[l][4].R0*Rz3(link[l][4].q);


        link[l][5].p = link[l][4].p + link[l][4].R*link[l][5].b;
//        printf("p5 = %.4f, %.4f, %.4f \n", link[l][5].p(0), link[l][5].p(1), link[l][5].p(2));
        link[l][5].R = link[l][4].R*link[l][5].R0*Rz3(link[l][5].q);

        link[l][6].p = link[l][5].p + link[l][5].R*link[l][6].b;
        link[l][6].R = link[l][5].R*link[l][6].R0*Rz3(link[l][6].q);

        vfo = link[l][6].R * Vec3<double>(1,0,0);
        //std::cout << "vfo = " << vfo << std::endl;
        tmp = vfo.transpose() * Vec3<double>(0,0,1);
        Terminal[l][4] = PI/2 - acos(tmp / vfo.norm());

        Terminal[l][0] = link[l][6].p(0);
        Terminal[l][1] = link[l][6].p(1);
        Terminal[l][2] = link[l][6].p(2);

        Terminal[l][3] = atan(sin(link[l][1].q)/cos(link[l][0].q)/cos(link[l][1].q));
        //Terminal[l].roll_ik = asin(sin(link[l][0].q)*cos(link[l][1].q));
    }
}
//data = x, y, z, yaw_ik, p, pas1, pas2
void LegClass::IKinematics(double* Tra, double* pas){
    int sign;
    double x, y, z, yaw_ik, p;
    double q[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Vec3<double> pos[7];
    Mat3<double> R[7];
    double tmp;

    for (int l = 0; l < 2; l ++) {
        sign =  (int)(1 - l*2);
        x = Tra[l*5 + 0];//-0.162944;
        y = Tra[l*5 + 1];//-0.104843;
        z = Tra[l*5 + 2];//-0.544781;
        yaw_ik = Tra[l*5 + 3];//0.356227;// = atan(sin(yaw_fk)/cos(roll_fk)/cos(yaw_fk));
        p = Tra[l*5 + 4];//deg2rad(10.0);// = p_fk;
        //q[4] = pas[l*2 + 0];//deg2rad(-3);
        //q[5] = pas[l*2 + 1];//deg2rad(5);

        pos[6] << x, y, z;

        pos[0] = link[l][0].b;
        R[0] = link[l][0].R0;
        //yaw_ik: #1世界坐标系 P~1, #0->#1
        Vec3<double> Ph1 = Rz3(-yaw_ik) * pos[6];

        // va = [y; z];
        // vb = [-z; y];
        // Pj1 = [ yj; zj];
        // Ph1 = [ yh; zh];
        // Pk1 = Pj1 + JK*va;
        // Pb1 = Pk1 + KB*vb;
        // bh1 = Ph1 - Pb1;
        // bh1y = bh1(1);%yh - yj - MK*y - JM*z + KB*z
        // bh1z = bh1(2);%zh - zj - MK*z + JM*y - KB*y
        // bh1//jk1
        // 1: bh1y/bh1z = (yh - yj - MK*y + (KB- JM)*z)/(zh - zj - MK*z - (KB - JM)*y) = y/z
        // 2: y^2 + z^2 = 1
        // right: 1->y = (yh - yj)/(zh - zj)*z + (KB- JM)/(zh - zj)  =kz + t
        // left : 1->y = (yh - yj)/(zh - zj)*z - (KB- JM)/(zh - zj)  =kz + t
        // yj = Pj1(2);
        // zj = Pj1(3);
        // yh = Ph1(2);
        // zh = Ph1(3);
        v1p2_linear eqt1;
        eqt1.h = 1;
//        binary_qua_cal.k = (Ph1(1) - Pj1(1)) / (Ph1(2) - Pj1(2));
//        binary_qua_cal.t = sign * (A.KB- A.JM) / (Ph1(2) - Pj1(2));
        eqt1.k = (Ph1(1) + A.LJ*sin(yaw_ik)) / (Ph1(2) - 0.0);
        eqt1.t = sign * (A.KB - A.MK) / (Ph1(2) - 0.0);
        eqt1.cal(1);

        double z1 = eqt1.output[0];
        double y1 = eqt1.output[1];
        //V2->V3 #1->#0
        Vec3<double> vb0 = Rz3(yaw_ik) * Vec3<double>(0, -z1, y1);
        Vec3<double> vb10 = Rz3(yaw_ik) * Vec3<double>(0, 1, 0);
        double tmp = vb0.transpose() * vb10;
        double roll_ik = acos(tmp / vb0.norm() / vb10.norm());
        if (y1 > 0) roll_ik = -roll_ik;
        q[0] = atan(sin(roll_ik) / cos(roll_ik) / cos(yaw_ik));
        q[1] = asin(sin(yaw_ik) * cos(roll_ik));

        R[0] = link[l][0].R0 * Rz3(q[0]);
        pos[1] = pos[0] + R[0] * link[l][1].b;//(A.OJ, 0, 0);
        R[1] = R[0] * link[l][1].R0 * Rz3(q[1]);
        pos[2] = pos[1] + R[1] * link[l][2].b;
        R[2] = R[1] * link[l][2].R0 * Rz3(q[2]);
        pos[3] = pos[2] + R[2] * link[l][3].b;
        //坐标系#2,->roll ->yaw,腿平面,x#3 = -x#0, y#3 = z#0
        Vec3<double> Pb2_3(0, 0, 0);
        Vec3<double> Ph2_3 = Rz3(-q[1]) * Rx3(q[0]) * (pos[6] - pos[3]);
        //v3->v2,
        Vec2<double> Pb2, Ph2;
        Pb2 << -Pb2_3(0), Pb2_3(2);
        Ph2 << -Ph2_3(0), Ph2_3(2);
        double h2 = abs(Ph2.norm());
        double theta_G = A.FGH - PI;//abs(q[5] + pas20) - abs(q[4] + pas10);
        double theta_D = A.CDG - PI;//abs(q[5] + pas20) - abs(q[4] + pas10);
        //double theta_D = (-1*sign*q[4] + sign*(link[l][3].q + sign*PI/2)) - PI;//A.FGH - PI;//abs(q[5] + pas20) - abs(q[4] + pas10);
        //double theta_G = sign*(q[5] + q[4]);//A.CDG - PI;//abs(q[5] + pas20) - abs(q[4] + pas10);
        //q[4] = -sign*(A.CDG - sign*q[3]);
        //q[5] = sign*(A.FGH - PI - sign*q[4]);
        //坐标系#3, pitch为0
        // 1 : (x + GH*sin(theta_G))^2 + (y + GH*cos(theta_G) + BD)^2 = DG^2
        // 2 : x^2 + y^2 = h2^2
        // 1 -> y = k*x + t
        v1p2_linear eqt2;
        eqt2.h = h2;
        eqt2.k = -2 * A.GH * sin(theta_G) / (2 * A.GH * cos(theta_G) + 2 * A.BD);
        eqt2.t = (A.DG * A.DG - h2 * h2 - A.GH * A.GH - A.BD * A.BD - 2 * A.BD * A.GH * cos(theta_G)) /
                           (2 * A.GH * cos(theta_G) + 2 * A.BD);
        eqt2.cal(0);

        double x3 = eqt2.output[0];
        double y3 = eqt2.output[1];

        Vec2<double> Ph3(x3, y3);
        Vec2<double> Pg3 = Ph3 + Rz(-theta_G) * Vec2<double>(0, A.GH);
        Vec2<double> Pb3(0, 0);
        Vec2<double> Pd3(0, -A.BD);
        Vec2<double> vgd3 = Pd3 - Pg3;
        Vec2<double> vba3 = Rz(-theta_D) * vgd3;
        Vec2<double> vbd3 = Pd3 - Pb3;
        tmp = vba3.transpose() * vbd3;
        q[3] = sign *(acos(tmp / vba3.norm() / vbd3.norm()) - PI/2);

        Vec2<double> vbh2 = Ph2 - Pb2;
        Vec2<double> vbh3 = Ph3 - Pb3;
        q[2] = sign * (asin(vbh3(0) / vbh3.norm()) - asin(vbh2(0) / vbh2.norm()));

        q[4] = -sign*(A.CDG - (sign*q[3]+PI/2));
        q[5] = sign*(A.FGH - PI - sign*q[4]);

        R[5] = link[l][0].R0*Rz3(q[0])*
               link[l][1].R0*Rz3(q[1])*
               link[l][2].R0*Rz3(q[2])*
               link[l][4].R0*Rz3(q[4])*
               link[l][5].R0*Rz3(q[5]);
        Mat3<double> tmpM;
        tmpM = R[5]*link[l][6].R0;
        v1p2 eqt3;
        double R31, R32, tmp_theta;
        R31 = tmpM(2, 0);
        R32 = tmpM(2, 1);
        tmp_theta = PI/2 - p;
        eqt3.a = R31*R31 + R32*R32;
        eqt3.b = -2*cos(tmp_theta)*R31;
        eqt3.c = cos(tmp_theta)*cos(tmp_theta) - R32*R32;
        //if (l == 0 & ((q[0] > 0 & q[1] > 0) | (q[0] < 0 & q[1] < 0)) |
        //l == 1 & ((q[0] < 0 & q[1] > 0) | (q[0] > 0 & q[1] < 0)))
        eqt3.cal(0);
        //else
        //eqt3.cal(1);
        q[6] = acos(eqt3.sol);

        double p_test, tmp_test;
        Vec3<double> vfo_test;
        vfo_test = R[5] * link[l][6].R0 * Rz3(q[6]) * Vec3<double>(1,0,0);
        tmp_test = vfo_test.transpose() * Vec3<double>(0,0,1);
        p_test = PI/2 - acos(tmp_test / vfo_test.norm());

        if (p_test - p > 0.00000000001 | p_test - p < -0.00000000001 )    q[6] = -q[6];

        for (int i = 0; i < 4; i ++)
            theta[l][i] = q[i];
        theta[l][4] = q[6];
        theta[l][5] = q[4];
        theta[l][6] = q[5];

    }
}




