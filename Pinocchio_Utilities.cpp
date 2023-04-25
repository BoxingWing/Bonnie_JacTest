//
// Created by boxingwang on 23-3-27.
//

#include "Pinocchio_Utilities.h"
const double pi=3.1415926;
Pinocchio_Utilities::Pinocchio_Utilities(std::string urdfName) {
    pinocchio::urdf::buildModel(urdfName,model_Bonnie_Static);
    pinocchio::JointModelFreeFlyer root_joint;
    pinocchio::urdf::buildModel(urdfName,root_joint,model_Bonnie_Dynamic);
    qB_urdf.setZero();
}

void Pinocchio_Utilities::setJointAngle(double *qr, double *ql, double *qPas_r, double *qPas_l) {
    qB_urdf<<ql[0],ql[1],ql[2],ql[3],qPas_l[0],qPas_l[1],ql[4],qr[0],qr[1],qr[2],qr[3],qPas_r[0],qPas_r[1],qr[4];
    qB_urdf(4)=qB_urdf(4)-98.66/180*pi;
    qB_urdf(5)=qB_urdf(5)-(-83.31/180*pi);
    qB_urdf(11)=qB_urdf(11)-(-98.66/180*pi);
    qB_urdf(12)=qB_urdf(12)-83.31/180*pi;
}

void Pinocchio_Utilities::computeJac() {

    pinocchio::Data data_B(model_Bonnie_Static);
    auto r_ankle_Joint=model_Bonnie_Static.getJointId("r_ankle_joint");
    auto l_ankle_Joint=model_Bonnie_Static.getJointId("l_ankle_joint");

    pinocchio::computeJointJacobians(model_Bonnie_Static,data_B,qB_urdf);
    Eigen::Matrix<double,4,7> J_L_tmp,J_R_tmp;
    Eigen::Matrix<double,6,14> J_tmp;
    J_tmp.setZero();
    pinocchio::getJointJacobian(model_Bonnie_Static,data_B,l_ankle_Joint,pinocchio::LOCAL_WORLD_ALIGNED,J_tmp);
    J_L_tmp.block<3,7>(0,0)=J_tmp.block<3,7>(0,0);
    J_L_tmp.row(3)=J_tmp.block<1,7>(5,0);
    J_tmp.setZero();
    pinocchio::getJointJacobian(model_Bonnie_Static,data_B,r_ankle_Joint,pinocchio::LOCAL_WORLD_ALIGNED,J_tmp);
    J_R_tmp.block<3,7>(0,0)=J_tmp.block<3,7>(0,7);
    J_R_tmp.row(3)=J_tmp.block<1,7>(5,7);

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

void Pinocchio_Utilities::computeIg() {
    pinocchio::Data data_B(model_Bonnie_Static);
    auto r_ankle_Joint=model_Bonnie_Static.getJointId("r_ankle_joint");
    auto l_ankle_Joint=model_Bonnie_Static.getJointId("l_ankle_joint");
    Eigen::VectorXd v_B = Eigen::VectorXd::Ones(model_Bonnie_Static.nv)*0;
    pinocchio::ccrba(model_Bonnie_Static,data_B,qB_urdf,v_B);
    Ig=data_B.Ig.inertia().matrix();
    pe_L=data_B.oMi[l_ankle_Joint].translation();
    pe_R=data_B.oMi[r_ankle_Joint].translation();
}

Eigen::Matrix<double, 3, 3> Pinocchio_Utilities::eul2Rot(double roll, double pitch, double yaw) {
    Eigen::Matrix<double,3,3> Rx,Ry,Rz;
    Rz<<cos(yaw),-sin(yaw),0,
        sin(yaw),cos(yaw),0,
        0,0,1;
    Ry<<cos(pitch),0,sin(pitch),
        0,1,0,
        -sin(pitch),0,cos(pitch);
    Rx<<1,0,0,
        0,cos(roll),-sin(roll),
        0,sin(roll),cos(roll);
    return Rz*Ry*Rx;
}

void Pinocchio_Utilities::computeG() {
    pinocchio::Data data(model_Bonnie_Static);
    pinocchio::computeGeneralizedGravity(model_Bonnie_Static,data,qB_urdf);
    Gq=data.g;
    //pinocchio::crba(model_Bonnie_Static,data,q_B);
}

double Pinocchio_Utilities::M8016_I2T(double Id)
{
    double Icmd=Id/66.0*4096.0;
    return Icmd*0.02245;
}

double Pinocchio_Utilities::M8016_T2I(double Td) {
    return Td/0.02245/4096.0*66.0;
}

double Pinocchio_Utilities::M10015_T2I(double Td) {
    double Icmd=(Td+ sgn(Td)*0.1141)/0.03189;
    return Icmd/4096.0*66.0;
}

double Pinocchio_Utilities::M10015_I2T(double Id)
{
    double Icmd=Id/66.0*4096.0;
    return Icmd*0.03189-sgn(Icmd)*0.1141;
}
double Pinocchio_Utilities::sgn(double in)
{
    if (in>0)
        return 1.0;
    else if (in<0)
        return -1.0;
    else
        return 0.0;
}