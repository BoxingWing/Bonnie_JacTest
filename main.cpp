#include <iostream>
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/jacobian.hpp"
#include "Leg.h"
#include "Pinocchio_Utilities.h"
#include "FileOperator.h"
#include "quill/Quill.h"

FileOperator fileRW("../TestData_Static_right.txt");
const double pi=3.1415926;
double Ts=0.001;
LegClass legKine;
Eigen::Matrix<double,6,1> FKDIY(Eigen::VectorXd q_B,bool removePasOff);
Eigen::Matrix<double,14,1> IKDIY(Eigen::Matrix<double,5,1> peR, Eigen::Matrix<double,5,1> peL,bool removePasOff);

int main() {
    Pinocchio_Utilities pinLib("../BonnieURDF_latest.urdf");
    const std::string urdf_filename_B = std::string("../BonnieURDF_latest.urdf");
    pinocchio::Model model_Bonnie;
    pinocchio::JointModelFreeFlyer root_joint;
//    pinocchio::urdf::buildModel(urdf_filename_B,root_joint,model_Bonnie);
    pinocchio::urdf::buildModel(urdf_filename_B,model_Bonnie);
    pinocchio::Data data_B(model_Bonnie);
    Eigen::VectorXd q_B = pinocchio::neutral(model_Bonnie);
    Eigen::VectorXd v_B = Eigen::VectorXd::Ones(model_Bonnie.nv)*0;
    Eigen::VectorXd a_B = Eigen::VectorXd::Ones(model_Bonnie.nv)*0;

    std::cout<<q_B.transpose()<<std::endl;

    double ql_robot[5],qr_robot[5],qPas_l[2],qPas_r[2];
    Eigen::Matrix<double,14,1> q_B_ori;
    q_B_ori<<-0.008,-0.0021,-0.5407,-0.4365,1.2890,-1.0082,0.1733,
            0.0033,0.0021,0.5582,0.4433,-1.2814,0.9995,-0.1779;
    qr_robot[0]=q_B_ori(0);qr_robot[1]=q_B_ori(1);qr_robot[2]=q_B_ori(2);qr_robot[3]=q_B_ori(3);qr_robot[4]=q_B_ori(6);
    qPas_r[0]=q_B_ori(4);qPas_r[5]=q_B_ori(5);
    ql_robot[0]=q_B_ori(7);ql_robot[1]=q_B_ori(8);ql_robot[2]=q_B_ori(9);ql_robot[3]=q_B_ori(10);ql_robot[4]=q_B_ori(13);
    qPas_l[0]=q_B_ori(11);qPas_r[0]=q_B_ori(12);
    q_B<<-0.008,-0.0021,-0.5407,-0.4365,1.2890,-1.0082,0.1733,
            0.0033,0.0021,0.5582,0.4433,-1.2814,0.9995,-0.1779;
    q_B(4)=q_B(4)-98.66/180*pi;
    q_B(5)=q_B(5)-(-83.31/180*pi);
    q_B(11)=q_B(11)-(-98.66/180*pi);
    q_B(12)=q_B(12)-83.31/180*pi;
    auto r_ankle_Joint=model_Bonnie.getJointId("r_ankle_joint");
    auto l_ankle_Joint=model_Bonnie.getJointId("l_ankle_joint");

    pinocchio::computeJointJacobians(model_Bonnie,data_B,q_B);
    Eigen::MatrixXd J_L(6,7);
    Eigen::MatrixXd J_tmp(6,14);
    J_tmp.setZero();
    pinocchio::getJointJacobian(model_Bonnie,data_B,l_ankle_Joint,pinocchio::LOCAL_WORLD_ALIGNED,J_tmp);
    J_L=J_tmp.block<6,7>(0,0);

    Eigen::Matrix<double,7,5> J_L_trans;
    Eigen::Matrix<double,6,5> J_L_active;
    J_L_trans.setZero();
    J_L_trans(0,0)=1;
    J_L_trans(1,1)=1;
    J_L_trans(2,2)=1;
    J_L_trans(3,3)=1;
    J_L_trans(4,3)=1;
    J_L_trans(5,3)=-1;
    J_L_trans(6,4)=1;

    J_L_active=J_L*J_L_trans;

    std::cout<<"J_L_trans"<<std::endl<<J_L_trans<<std::endl;
    std::cout<<"J_L"<<std::endl<<J_L<<std::endl;
    std::cout<<"J_L_active"<<std::endl<<J_L_active<<std::endl;
    std::cout<<data_B.oMi[l_ankle_Joint].translation().transpose()<<std::endl;

    Eigen::Matrix<double,6,1> FeDIY_nominal,FeDIY_tmp;
    FeDIY_nominal=FKDIY(q_B,true);
    //std::cout<<FeDIY.transpose()<<std::endl;
    Eigen::Matrix<double,3,5> J_L_DIY;
    Eigen::Matrix<double,3,1> J_L_col;
    double deltaAngle=0.001;
    for (int i=0;i<5;i++){
        auto q_tmp=q_B;
        if (i<4)
            q_tmp(i)+=deltaAngle;
        else
            q_tmp(6)+=deltaAngle;
        FeDIY_tmp= FKDIY(q_tmp,true);
        J_L_col<<FeDIY_tmp(3)-FeDIY_nominal(3),FeDIY_tmp(4)-FeDIY_nominal(4),FeDIY_tmp(5)-FeDIY_nominal(5);
        J_L_col=J_L_col/deltaAngle;
        J_L_DIY.col(i)=J_L_col;
    }
    std::cout<<"J_L_DIY"<<std::endl;
    std::cout<<J_L_DIY<<std::endl;
    // test for the change of Ig under different leg length
    pinLib.setJointAngle(qr_robot,ql_robot,qPas_r,qPas_l);
    pinLib.computeJac();
    std::cout<<pinLib.J_L<<std::endl;
    std::cout<<pinLib.J_R<<std::endl;

    pinLib.computeIg();
    std::cout<<pinLib.Ig<<std::endl;

    //================ Test for jacobian estimated motor torques ==========
//    std::vector<double> tmpValue;
//    std::string tmpStr;
//    quill::Handler *file_handler = quill::file_handler("../Outputdata_right.txt", "w");  // Get a pointer to the default logger
//    file_handler->set_pattern(QUILL_STRING("%(message)"));
//    quill::Logger *dl = quill::create_logger("logger", file_handler);
//    quill::start();
//    double ql[5],qr[5],pas[4],Il[4],Ir[5];
//    double acc[3],eul[3],omegaL[3],phaseAll,legSwingInd,pB_fk[6],pas_delta[4];
//    Eigen::Matrix<double,3,3> Reul;
//    for (int i = 0; i < fileRW.getTotalLine(); i++) {
//        // read data
//        fileRW.getNewLine();
//        fileRW.getNumsInLine();
//        acc[0] = fileRW.values[0];
//        acc[1] = fileRW.values[1];
//        acc[2] = fileRW.values[2];
//        eul[0] = fileRW.values[3];
//        eul[1] = fileRW.values[4];
//        eul[2] = fileRW.values[5];
//        omegaL[0] = fileRW.values[6];
//        omegaL[1] = fileRW.values[7];
//        omegaL[2] = fileRW.values[8];
//        phaseAll = fileRW.values[9];
//        legSwingInd = fileRW.values[10];
//        pB_fk[0] = fileRW.values[11];
//        pB_fk[1] = fileRW.values[12];
//        pB_fk[2] = fileRW.values[13];
//        pB_fk[3] = fileRW.values[14];
//        pB_fk[4] = fileRW.values[15];
//        pB_fk[5] = fileRW.values[16];
//        pas_delta[0] = fileRW.values[17];
//        pas_delta[1] = fileRW.values[18];
//        pas_delta[2] = fileRW.values[19];
//        pas_delta[3] = fileRW.values[20];
//        qr[0]=fileRW.values[21];qr[1]=fileRW.values[22];qr[2]=fileRW.values[23];qr[3]=fileRW.values[24];qr[4]=fileRW.values[25];
//        pas[0]=fileRW.values[26];pas[1]=fileRW.values[27];
//        ql[0]=fileRW.values[28];ql[1]=fileRW.values[29];ql[2]=fileRW.values[30];ql[3]=fileRW.values[31];ql[4]=fileRW.values[32];
//        pas[2]=fileRW.values[33];pas[3]=fileRW.values[34];
//        Ir[0]=fileRW.values[35];Ir[1]=fileRW.values[36];Ir[2]=fileRW.values[37];Ir[3]=fileRW.values[38];
//        Il[0]=fileRW.values[39];Il[1]=fileRW.values[40];Il[2]=fileRW.values[41];Il[3]=fileRW.values[42];
//
//        q_B_ori<<ql[0],ql[1],ql[2],ql[3],pas[2],pas[3],ql[4],qr[0],qr[1],qr[2],qr[3],pas[0],pas[1],qr[4];
//        qPas_r[0]=pas[0];qPas_r[1]=pas[1];qPas_l[0]=pas[2];qPas_l[1]=pas[3];
//        pinLib.setJointAngle(qr,ql,qPas_r,qPas_l);
//        pinLib.computeJac();
//        //std::cout<<pinLib.J_L<<std::endl;
//        //std::cout<<pinLib.J_R<<std::endl;
//        Reul=Pinocchio_Utilities::eul2Rot(eul[0],eul[1],0);
//        Eigen::Matrix<double,5,1> tauL,tauR,tauL_M,tauR_M;
//        Eigen::Matrix<double,3,1> FendL,FendR;
//        FendL<<-0.4212, -3.8698, -70.7308;
//        FendR<<-3.1929, 5.9263, -54.5708;
//
//        Eigen::Matrix<double,4,1> WrenchL,WrenchR;
//        Eigen::Matrix<double,3,1> FendLW,FendRW;
//        FendLW=Reul.transpose()*FendL;
//        FendRW=Reul.transpose()*FendR;
//
//        WrenchL<<FendLW(0),FendLW(1),FendLW(2),0;
//        WrenchR<<FendRW(0),FendRW(1),FendRW(2),0;
//
//        tauL=pinLib.J_L.transpose()*WrenchL;
//        tauR=pinLib.J_R.transpose()*WrenchR;
//        tauL_M<< pinLib.M10015_I2T(Il[0]),pinLib.M8016_I2T(Il[1]),pinLib.M8016_I2T(Il[2]),pinLib.M10015_I2T(Il[3]),0;
//        tauR_M<< pinLib.M10015_I2T(Ir[0]),pinLib.M8016_I2T(Ir[1]),pinLib.M8016_I2T(Ir[2]),pinLib.M10015_I2T(Ir[3]),0;
//
//        std::cout<<"================"<<std::endl;
//        std::cout<<"pitch angle: "<<eul[1]<<std::endl;
//        std::cout<<"LLeg_Jac:"<<tauL.transpose()<<std::endl;
//        std::cout<<"LLeg_Mot:"<<tauL_M.transpose()<<std::endl<<std::endl;
//        std::cout<<"RLeg_Jac:"<<tauR.transpose()<<std::endl;
//        std::cout<<"RLeg_Mot:"<<tauR_M.transpose()<<std::endl;
//
//        tmpValue.clear();
//
//        for (int ii=0;ii<5;ii++)
//            tmpValue.push_back(tauL(ii));
//        for (int ii=0;ii<5;ii++)
//            tmpValue.push_back(tauL_M(ii));
//        for (int ii=0;ii<5;ii++)
//            tmpValue.push_back(tauR(ii));
//        for (int ii=0;ii<5;ii++)
//            tmpValue.push_back(tauR_M(ii));
//
//        tmpStr = fmt::format("{:.5f}", fmt::join(tmpValue, ","));
//        LOG_INFO(dl, "{}", tmpStr);
//    }
//    std::cout<<std::endl;
//    std::cout<<pinLib.pe_L.transpose()<<std::endl;
//    std::cout<<pinLib.pe_R.transpose()<<std::endl;
//
//    pinLib.computeG();
//    std::cout<<pinLib.Gq.transpose()<<std::endl;

    //---------------------- Test for centroid inertia -------------------------
    Eigen::Matrix<double,14,1> qIK;
    Eigen::Matrix<double,5,1> peL,peR;
    peL<<0,0.0824,-0.6,0,0;
    peR<<0,-0.0824,-0.6,0,0;
    qIK=IKDIY( peR, peL, true);
//    Eigen::Matrix<double,6,1> peFK;
//    peFK=FKDIY(qIK);
//    std::cout<<peFK.transpose()<<std::endl;
    std::cout<<qIK.transpose()<<std::endl;

    pinLib.qB_urdf=qIK;
    pinLib.computeIg();
    std::cout<<"Ig"<<std::endl;
    std::cout<<pinLib.Ig<<std::endl;
    std::cout<<"pe_R "<<pinLib.pe_R.transpose()<<std::endl;
    std::cout<<"pe_L "<<pinLib.pe_L.transpose()<<std::endl;
    std::cout<<"pCoM "<<pinLib.pCoM.transpose()<<std::endl;
    std::cout<<"totalMass "<<pinLib.totalMass<<std::endl;
    pinLib.computeG();
    std::cout<<"Gq "<<pinLib.Gq.transpose()<<std::endl;

}

Eigen::Matrix<double,6,1> FKDIY(Eigen::VectorXd q_B,bool removePasOff)
{
    auto q_kine=q_B;
    double tmp[3];
    tmp[0]=q_kine(4);tmp[1]=q_kine(5);tmp[2]=q_kine(6);
    if (removePasOff)
    {   tmp[0]=tmp[0]+(98.66/180*pi);
        tmp[1]=tmp[1]+(-83.31/180*pi);}
    q_kine(4)=tmp[2];q_kine(5)=tmp[0];q_kine(6)=tmp[1];

    tmp[0]=q_kine(11);tmp[1]=q_kine(12);tmp[2]=q_kine(13);
    if (removePasOff) {
        tmp[0] = tmp[0] + (-98.66 / 180 * pi);
        tmp[1] = tmp[1] + (83.31 / 180 * pi);
    }
    q_kine(11)=tmp[2];q_kine(12)=tmp[0];q_kine(13)=tmp[1];
    for (int i=0;i<7;i++)
    {
        double tmp2=0;
        tmp2=q_kine(i);
        q_kine(i)=q_kine(i+7);
        q_kine(i+7)=tmp2;
    }
    legKine.FKinematics(q_kine.data());
    //std::cout<<legKine.Terminal[1][0]<<' '<<legKine.Terminal[1][1]+0.125<<' '<<legKine.Terminal[1][2];
    Eigen::Matrix<double,6,1> Fe;
    Fe<<legKine.Terminal[0][0],legKine.Terminal[0][1]-0.125,legKine.Terminal[0][2],
        legKine.Terminal[1][0],legKine.Terminal[1][1]+0.125,legKine.Terminal[1][2];
    return Fe;
}

Eigen::Matrix<double,14,1> IKDIY(Eigen::Matrix<double,5,1> peR, Eigen::Matrix<double,5,1> peL, bool removePasOff)
{
    double pas[4]={0};
    double pe[10]={0};
    peR(1)=peR(1)+0.125;
    peL(1)=peL(1)-0.125;
    for (int i=0;i<5;i++)
        {pe[i]=peR(i);pe[i+5]=peL(i);}
    legKine.IKinematics(pe,pas);
    Eigen::Matrix<double,14,1> qFK;
//    qFK = ql[0],ql[1],ql[2],ql[3],pas[2],pas[3],ql[4],qr[0],qr[1],qr[2],qr[3],pas[0],pas[1],qr[4];
    qFK<<legKine.theta[1][0],legKine.theta[1][1],legKine.theta[1][2],legKine.theta[1][3],legKine.theta[1][5],legKine.theta[1][6], legKine.theta[1][4],
            legKine.theta[0][0],legKine.theta[0][1],legKine.theta[0][2],legKine.theta[0][3],legKine.theta[0][5],legKine.theta[0][6], legKine.theta[0][4];
    if (removePasOff)
    {
        qFK(4)=qFK(4)-(98.66/180*pi);
        qFK(5)=qFK(5)-(-83.31/180*pi);
        qFK(11)=qFK(11)-(-98.66/180*pi);
        qFK(12)=qFK(12)-(83.31/180*pi);
    }

    return qFK;
}