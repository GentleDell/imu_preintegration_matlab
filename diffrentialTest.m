clear;clc;
% global attiCaculator;

randn('seed',0);
step = 0.01;
start_time = 0;
end_time = 50;
tspan = [start_time:step:end_time]';
N = length(tspan);
Ar = 10;
r = [Ar*sin(tspan) Ar*cos(tspan) 0.5*tspan.*tspan];
v = [Ar*cos(tspan) -Ar*sin(tspan) tspan];
acc_inertial = [-Ar*sin(tspan) -Ar*cos(tspan) ones(N,1)];

atti = [0.1*sin(tspan) 0.1*sin(tspan) 0.1*sin(tspan)];
Datti = [0.1*cos(tspan) 0.1*cos(tspan) 0.1*cos(tspan)];
g = [0 0 -9.8]';
gyro_pure = zeros(N,3);
acc_pure = zeros(N,3);

a = wgn(N,1,1)/5;
b = zeros(N,1);
b(1) = a(1)*step;
gt = zeros(3,4,N);
%生成无噪声的imu数据
for iter = 1:N
    A = AttitudeBase.Datti2w(atti(iter,:));
    gyro_pure(iter,:) = Datti(iter,:)*A';
    cnb = AttitudeBase.a2cnb(atti(iter,:));
    acc_pure(iter,:) = cnb*(acc_inertial(iter,:)' - g);
    gt(:,:,iter) = [cnb' r(iter,:)'];
end

[acc2,gyro2] = IMURawData.DifferentialIMUKinectics(gt,step);