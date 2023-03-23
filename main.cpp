#include <iostream>
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "Leg.h"
const double pi=3.1415926;
double Ts=0.001;
LegClass legKine;

int main() {
    const std::string urdf_filename_B = std::string("../BonnieURDF_latest.urdf");
    pinocchio::Model model_Bonnie;
    pinocchio::JointModelFreeFlyer root_joint;
    //pinocchio::urdf::buildModel(urdf_filename_B,root_joint,model_Bonnie);
    pinocchio::urdf::buildModel(urdf_filename_B,model_Bonnie);
    pinocchio::Data data_B(model_Bonnie);
    Eigen::VectorXd q_B = pinocchio::neutral(model_Bonnie);
    Eigen::VectorXd v_B = Eigen::VectorXd::Ones(model_Bonnie.nv)*0;
    Eigen::VectorXd a_B = Eigen::VectorXd::Ones(model_Bonnie.nv)*0;
    /*q_B<<0,0,0,
            0,0,0,1,
            -0.008,-0.0021,-0.5407,-0.4365,1.2890,-1.0082,0.1733,
            0.0033,0.0021,0.5582,0.4433,-1.2814,0.9995,-0.1779;
    q_B(11)=q_B(11)-98.66/180*pi;
    q_B(12)=q_B(12)-(-83.31/180*pi);
    q_B(18)=q_B(18)-(-98.66/180*pi);
    q_B(19)=q_B(19)-83.31/180*pi;*/
    q_B<<-0.008,-0.0021,-0.5407,-0.4365,1.2890,-1.0082,0.1733,
            0.0033,0.0021,0.5582,0.4433,-1.2814,0.9995,-0.1779;
    q_B(4)=q_B(4)-98.66/180*pi;
    q_B(5)=q_B(5)-(-83.31/180*pi);
    q_B(11)=q_B(11)-(-98.66/180*pi);
    q_B(12)=q_B(12)-83.31/180*pi;
    auto r_ankle_Joint=model_Bonnie.getJointId("r_ankle_joint");
    auto l_ankle_Joint=model_Bonnie.getJointId("l_ankle_joint");

    pinocchio::computeJointJacobians(model_Bonnie,data_B,q_B);
    Eigen::MatrixXd J_L(6,model_Bonnie.nv);
    J_L.setZero();
    pinocchio::getJointJacobian(model_Bonnie,data_B,l_ankle_Joint,pinocchio::LOCAL_WORLD_ALIGNED,J_L);
    std::cout<<J_L<<std::endl;
    std::cout<<J_L.cols()<<std::endl;
    std::cout<<data_B.oMi[l_ankle_Joint].translation().transpose()<<std::endl;
    //std::cout<<data_B.J<<std::endl;

    auto q_kine=q_B;
    double tmp[3];
    tmp[0]=q_kine(4);tmp[1]=q_kine(5);tmp[2]=q_kine(6);
    tmp[0]=tmp[0]+(98.66/180*pi);
    tmp[1]=tmp[1]+(-83.31/180*pi);
    q_kine(4)=tmp[2];q_kine(5)=tmp[0];q_kine(6)=tmp[1];

    tmp[0]=q_kine(11);tmp[1]=q_kine(12);tmp[2]=q_kine(13);
    tmp[0]=tmp[0]+(-98.66/180*pi);
    tmp[1]=tmp[1]+(83.31/180*pi);
    q_kine(11)=tmp[2];q_kine(12)=tmp[0];q_kine(13)=tmp[1];
    for (int i=0;i<7;i++)
    {
        double tmp2=0;
        tmp2=q_kine(i);
        q_kine(i)=q_kine(i+7);
        q_kine(i+7)=tmp2;
    }
    legKine.FKinematics(q_kine.data());
    std::cout<<legKine.Terminal[1][0]<<' '<<legKine.Terminal[1][1]+0.125<<' '<<legKine.Terminal[1][2];
}
