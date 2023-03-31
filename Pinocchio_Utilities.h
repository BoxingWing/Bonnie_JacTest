//
// Created by boxingwang on 23-3-27.
//

#ifndef BONNIE_JACTEST_PINOCCHIO_UTILITIES_H
#define BONNIE_JACTEST_PINOCCHIO_UTILITIES_H

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "pinocchio/algorithm/centroidal.hpp"

class Pinocchio_Utilities {
public:
    pinocchio::Model model_Bonnie_Static, model_Bonnie_Dynamic;
    Eigen::VectorXd q_S, q_D;
    Eigen::VectorXd dq_S, dq_D;
    Eigen::Matrix<double,3,5> J_L,J_R;
    Eigen::Matrix<double,3,3> Ig;
    Eigen::Matrix<double,3,1> pe_L,pe_R;

    Pinocchio_Utilities(std::string urdfName);

    // NOTE: q_B here should be the joint angles whose offset are defined in the real system, NOT in the urdf
    // HOWEVER, the joint order in q_B should follow the one defined in the urdf
    void computeIg(Eigen::VectorXd q_B);
    void computeJac(Eigen::VectorXd q_B);
};


#endif //BONNIE_JACTEST_PINOCCHIO_UTILITIES_H
