//
// Created by boxingwang on 23-3-27.
//

#include "Pinocchio_Utilities.h"
const double pi=3.1415926;
Pinocchio_Utilities::Pinocchio_Utilities(std::string urdfName) {
    pinocchio::urdf::buildModel(urdfName,model_Bonnie_Static);
    pinocchio::JointModelFreeFlyer root_joint;
    pinocchio::urdf::buildModel(urdfName,root_joint,model_Bonnie_Dynamic);
}


void Pinocchio_Utilities::computeJac(Eigen::VectorXd q_B) {
    pinocchio::Data data_B(model_Bonnie_Static);
    q_B(4)=q_B(4)-98.66/180*pi;
    q_B(5)=q_B(5)-(-83.31/180*pi);
    q_B(11)=q_B(11)-(-98.66/180*pi);
    q_B(12)=q_B(12)-83.31/180*pi;
    auto r_ankle_Joint=model_Bonnie_Static.getJointId("r_ankle_joint");
    auto l_ankle_Joint=model_Bonnie_Static.getJointId("l_ankle_joint");

    pinocchio::computeJointJacobians(model_Bonnie_Static,data_B,q_B);
    Eigen::Matrix<double,3,7> J_L_tmp,J_R_tmp;
    Eigen::Matrix<double,6,14> J_tmp;
    J_tmp.setZero();
    pinocchio::getJointJacobian(model_Bonnie_Static,data_B,l_ankle_Joint,pinocchio::LOCAL_WORLD_ALIGNED,J_tmp);
    J_L_tmp=J_tmp.block<3,7>(0,0);
    J_tmp.setZero();
    pinocchio::getJointJacobian(model_Bonnie_Static,data_B,r_ankle_Joint,pinocchio::LOCAL_WORLD_ALIGNED,J_tmp);
    J_R_tmp=J_tmp.block<3,7>(0,7);

    Eigen::Matrix<double,7,5> J_trans;
    J_trans.setZero();
    J_trans(0,0)=1;
    J_trans(1,1)=1;
    J_trans(2,2)=1;
    J_trans(3,3)=1;
    J_trans(4,3)=1;
    J_trans(5,3)=-1;
    J_trans(6,4)=1;

    J_L=J_L_tmp*J_trans;
    J_R=J_R_tmp*J_trans;
    pe_L=data_B.oMi[l_ankle_Joint].translation();
    pe_R=data_B.oMi[r_ankle_Joint].translation();
}

void Pinocchio_Utilities::computeIg(Eigen::VectorXd q_B) {
    pinocchio::Data data_B(model_Bonnie_Static);
    Eigen::VectorXd v_B = Eigen::VectorXd::Ones(model_Bonnie_Static.nv)*0;
    q_B(4)=q_B(4)-98.66/180*pi;
    q_B(5)=q_B(5)-(-83.31/180*pi);
    q_B(11)=q_B(11)-(-98.66/180*pi);
    q_B(12)=q_B(12)-83.31/180*pi;
    pinocchio::ccrba(model_Bonnie_Static,data_B,q_B,v_B);
    Ig=data_B.Ig.inertia().matrix();
}