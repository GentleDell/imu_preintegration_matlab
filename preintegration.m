clear;clc;
% global attiCaculator;

randn('seed',0);
step = 0.01;
start_time = 0;
end_time = 50;
tspan = [start_time:step:end_time]';
N = length(tspan);
Ar = 10;
r = [Ar*sin(tspan) Ar*cos(tspan) 0.5*tspan.*tspan]; % pose of GPS —— global pose
v = [Ar*cos(tspan) -Ar*sin(tspan) tspan];
acc_inertial = [-Ar*sin(tspan) -Ar*cos(tspan) ones(N,1)];  % global acceleration

atti = [0.1*sin(tspan) 0.1*sin(tspan) 0.1*sin(tspan)];  % Rotation: pitch roll yaw || x y z
Datti = [0.1*cos(tspan) 0.1*cos(tspan) 0.1*cos(tspan)];
g = [0 0 -9.8]';
gyro_pure = zeros(N,3);
acc_pure = zeros(N,3);

a = wgn(N,1,1)/5;
b = zeros(N,1);
b(1) = a(1)*step;
%% generate imu data
for iter = 1:N
    A = AttitudeBase.Datti2w(atti(iter,:));
    gyro_pure(iter,:) = Datti(iter,:)*A';
    cnb = AttitudeBase.a2cnb(atti(iter,:));
    acc_pure(iter,:) = cnb*(acc_inertial(iter,:)' - g);
end
% state0 = zeros(10,1);  
% state0(7) = 1;

% add noise and bias (0.5 m/s^2 & 1 degree)
acc_noised = acc_pure + 0.1*randn(N,3) + 0.5*ones(N,3);
gyro_noised = gyro_pure + (randn(N,3)*2 + ones(N,3))/180*pi;
accCov = 0.01*eye(3);
gyroCov = (2/180*pi)^2*eye(3);

%% preintegration
imuPara = IMUPara(accCov,gyroCov,ones(3),ones(3));  % initial imuPara with ones?

PIM = PreintegrateMeasurement();
for i = 1:100   % the former 100 imu data
    PIM = PIM.Preintegrate(acc_noised(i,:)',gyro_noised(i,:)',imuPara,step);
end

si =  NavState(zeros(3,1),r(1,:)',v(1,:)',zeros(3,1),zeros(3,1));   % initial the NavState
sj = si.predict(PIM);   % update the NavState by imu preintegraion

% sj = si;
% cov = zeros(15,15);
% cov(1:9,1:9) = PIM.cov_;
% cov(10:15,10:15) = 0.1*eye(6);
% infor = eye(15)/cov;
% oriMeas = SO3.log(AttitudeBase.a2cnb(atti(100,:))');

cov = zeros(12,12);
cov(1:9,1:9) = PIM.cov_;
cov(10:12,10:12) = 0.1*eye(3);
infor = eye(12)/cov;
for i = 1:2
    %     [ err,J ] = IMUErrorJacobian.GPSAndOri( si,sj,PIM,r(100,:)',oriMeas );
    [ err,J ] = IMUErrorJacobian.GPS( si,sj,PIM,r(100,:)' );    %% get residual error & uncomplete Jaccobian
    fprintf('err: %f\n',norm(err));
    %     if norm(err)<0.00001
    %         break;
    %     end
    delta = -(J'*infor*J)\J'*infor*err;
    %     fprintf('delta: %f\n',norm(delta))
    
    sj = sj.update(delta);
    si.ba_ = sj.ba_;
    si.bg_ = sj.bg_;
end



