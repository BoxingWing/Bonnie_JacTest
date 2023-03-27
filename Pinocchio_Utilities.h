//
// Created by boxingwang on 23-3-27.
//

#ifndef BONNIE_JACTEST_PINOCCHIO_UTILITIES_H
#define BONNIE_JACTEST_PINOCCHIO_UTILITIES_H

#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"

class Pinocchio_Utilities {
public:
    pinocchio::Model model_Bonnie_Static, model_Bonnie_Dynamic;
    Eigen::VectorXd q_S, q_D;
    Eigen::VectorXd dq_S, dq_D;
    Eigen::Matrix<double,3,5> J_L,J_R;

    Pinocchio_Utilities(std::string urdfName);
    void computeJac(Eigen::VectorXd q_B);
};


#endif //BONNIE_JACTEST_PINOCCHIO_UTILITIES_H
