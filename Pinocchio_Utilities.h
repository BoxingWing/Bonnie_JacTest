//
// Created by boxingwang on 23-3-27.
//

#ifndef BONNIE_JACTEST_PINOCCHIO_UTILITIES_H
#define BONNIE_JACTEST_PINOCCHIO_UTILITIES_H

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/rnea.hpp"

class Pinocchio_Utilities {
public:
    pinocchio::Model model_Bonnie_Static, model_Bonnie_Dynamic;
    Eigen::VectorXd q_S, q_D;
    Eigen::VectorXd dq_S, dq_D;
    Eigen::Matrix<double,4,5> J_L,J_R;
    Eigen::Matrix<double,3,3> Ig;
    Eigen::Matrix<double,3,1> pe_L,pe_R;
    Eigen::Matrix<double,14,1> Gq; // generalized gravity term
    Eigen::Matrix<double,14,1> qB_urdf;

    Pinocchio_Utilities(std::string urdfName);

    // NOTE: qr, ql, qPas_r, qPas_l should be direct values get from the robot
    void setJointAngle(double* qr, double *ql, double *qPas_r, double *qPas_l);
    void computeIg();
    void computeJac();
    void computeG();
    double M8016_I2T(double Id);
    double M10015_I2T(double Id);
    double M8016_T2I(double Td);
    double M10015_T2I(double Td);
    double sgn(double in);

    static Eigen::Matrix<double,3,3> eul2Rot(double roll, double pitch, double yaw);
};


#endif //BONNIE_JACTEST_PINOCCHIO_UTILITIES_H
