clear variables;
close all;
data=load('quill_outputdata.txt');
tauL=data(:,1:5);
tauL_M=data(:,6:10);
tauR=data(:,11:15);
tauR_M=data(:,16:20);
figure('name','Left leg');
subplot(2,3,1)
plot(tauL(:,1));hold on;
plot(tauL_M(:,1));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauL_M(:,1)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,2)
plot(tauL(:,2));hold on;
plot(tauL_M(:,2));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauL_M(:,2)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,3)
plot(tauL(:,3));hold on;
plot(tauL_M(:,3));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauL_M(:,3)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava')
subplot(2,3,4)
plot(tauL(:,4));hold on;
plot(tauL_M(:,4));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauL_M(:,4)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,5)
plot(tauL(:,5));hold on;
plot(tauL_M(:,5));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauL_M(:,5)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');

figure('name','Right leg');
subplot(2,3,1)
plot(tauR(:,1));hold on;
plot(tauR_M(:,1));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauR_M(:,1)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,2)
plot(tauR(:,2));hold on;
plot(tauR_M(:,2));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauR_M(:,2)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,3)
plot(tauR(:,3));hold on;
plot(tauR_M(:,3));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauR_M(:,3)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,4)
plot(tauR(:,4));hold on;
plot(tauR_M(:,4));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauR_M(:,4)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');
subplot(2,3,5)
plot(tauR(:,5));hold on;
plot(tauR_M(:,5));
plot([1,length(tauL_M(:,1))],[1,1]*mean(tauR_M(:,5)));
legend('Jac\_Est','Mot\_Read','Mot\_Read\_ava');






