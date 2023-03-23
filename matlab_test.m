clear variables;
close all;
b=importrobot('BonnieURDF_latest.urdf');
config = homeConfiguration(b);
config(1).JointPosition=-0.008;
config(2).JointPosition=-0.0021;
config(3).JointPosition=-0.5407;
config(4).JointPosition=-0.4365;
config(5).JointPosition=1.2890-98.66/180*pi;
config(6).JointPosition=-1.0082-(-83.31/180*pi);
config(7).JointPosition=0.1733;
J=geometricJacobian(b,config,'l_ankle_link')
show(b,config)
res=getTransform(b,config,'l_ankle_link')

