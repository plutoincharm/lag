clearvars; clc;

global m L r MI T2a tf

L(1) = 3;
L(2) = 7;
L(3) = 6;
r(1) = L(1)/2;
r(2) = L(2)/2;
r(3) = L(3)/2;
m(1) = 1;
m(2) = 3;
m(3) = 2;
MI(1) = 1/12*m(1)*L(1)^2;
MI(2) = 1/12*m(2)*L(2)^2;
MI(3) = 1/12*m(3)*L(3)^2;
T2a = 2;
th2 = pi/3;
th3 = 0.380251206692934;%0.3803;
th4 = pi/3;

tf = 12;
dt = 1e-3;
t = 0:dt:tf;
n = length(t);
%%%%%%% var order: (omg2,omg3,omg4,th2,th3,th4)
q0 = [0;0;0;th2;th3;th4];

handle_f = @fun4barLagcODE;

%%%%%%%%%%% Explicit Euler %%%%%%%%%%%%%%%%%
for i = 1:n-1
    if i==1
        q(i,:) = q0';
        qnew = fun_explicitEuler(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    else
        qnew = fun_explicitEuler(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    end
end
qExEu = q;
% save('fourbarNwtODE_ExEu.mat','t','q');
omg2ExEu = qExEu(:,1);
omg3ExEu = qExEu(:,2);
omg4ExEu = qExEu(:,3);
tht2ExEu = qExEu(:,4);
tht3ExEu = qExEu(:,5);
tht4ExEu = qExEu(:,6);

figure; hold on; grid on;
plot(t,tht2ExEu);
plot(t,tht3ExEu);
plot(t,tht4ExEu);
legend('tht2','tht3','tht4');
title('Explicit Euler');

figure; hold on; grid on;
plot(t,omg2ExEu);
plot(t,omg3ExEu);
plot(t,omg4ExEu);
legend('omg2','omg3','omg4');
title('Explicit Euler');

%%%%%%%%%%% Predictor-Corrector %%%%%%%%%%%%%%%%%
for i = 1:n-1
    if i==1
        q(i,:) = q0';
        qnew = fun_PredCor(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    else
        qnew = fun_PredCor(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    end
end
qPC = q;
% save('fourbarNwtODE_PC.mat','t','q');
omg2PC = qPC(:,1);
omg3PC = qPC(:,2);
omg4PC = qPC(:,3);
tht2PC = qPC(:,4);
tht3PC = qPC(:,5);
tht4PC = qPC(:,6);

figure; hold on; grid on;
plot(t,tht2PC);
plot(t,tht3PC);
plot(t,tht4PC);
legend('tht2','tht3','tht4');
title('Predictor Corrector');

figure; hold on; grid on;
plot(t,omg2PC);
plot(t,omg3PC);
plot(t,omg4PC);
legend('omg2','omg3','omg4');
title('Predictor Corrector');

%%%%%%%%%%% Huen %%%%%%%%%%%%%%%%%
for i = 1:n-1
    if i==1
        q(i,:) = q0';
        qnew = fun_Huen(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    else
        qnew = fun_Huen(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    end
end
qHuen = q;
% save('fourbarNwtODE_Huen.mat','t','q');
omg2Huen = qHuen(:,1);
omg3Huen = qHuen(:,2);
omg4Huen = qHuen(:,3);
tht2Huen = qHuen(:,4);
tht3Huen = qHuen(:,5);
tht4Huen = qHuen(:,6);

figure; hold on; grid on;
plot(t,tht2Huen);
plot(t,tht3Huen);
plot(t,tht4Huen);
legend('tht2','tht3','tht4');
title('Huen');

figure; hold on; grid on;
plot(t,omg2Huen);
plot(t,omg3Huen);
plot(t,omg4Huen);
legend('omg2','omg3','omg4');
title('Huen');

%%%%%%%%%%% RK-4 %%%%%%%%%%%%%%%%%
for i = 1:n-1
    if i==1
        q(i,:) = q0';
        qnew = fun_rk4(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    else
        qnew = fun_rk4(handle_f,q(i,:)',dt,t(i));
        q(i+1,:) = qnew';
    end
end
qRK4 = q;
% save('fourbarNwtODE_RK4.mat','t','q');
omg2RK4 = qRK4(:,1);
omg3RK4 = qRK4(:,2);
omg4RK4 = qRK4(:,3);
tht2RK4 = qRK4(:,4);
tht3RK4 = qRK4(:,5);
tht4RK4 = qRK4(:,6);

figure; hold on; grid on;
plot(t,tht2RK4);
plot(t,tht3RK4);
plot(t,tht4RK4);
legend('tht2','tht3','tht4');
title('RK-4');

figure; hold on; grid on;
plot(t,omg2RK4);
plot(t,omg3RK4);
plot(t,omg4RK4);
legend('omg2','omg3','omg4');
title('RK-4');
