

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
%生成无噪声的imu数据
for iter = 1:N
    A = AttitudeBase.Datti2w(atti(iter,:));
    gyro_pure(iter,:) = Datti(iter,:)*A';
    cnb = AttitudeBase.a2cnb(atti(iter,:));
    acc_pure(iter,:) = cnb*(acc_inertial(iter,:)' - g);
end
% state0 = zeros(10,1);
% state0(7) = 1;

%加速度计和陀螺仪加噪声
acc_noise = acc_pure + 0.1*randn(N,3) + 0.5*ones(N,3);
gyro_noise = gyro_pure + (randn(N,3)*2 + ones(N,3))/180*pi;

accCov = 0.01*ones(3);
gyroCov = (2/180*pi)^2*ones(3);

imuPara = IMUPara(accCov,gyroCov,ones(3),ones(3));

PIM = PreintegrateMeasurement();
for i = 1:100
    PIM = PIM.Preintegrate(acc_noise(i,:)',gyro_noise(i,:)',imuPara,step);
end

pErrors_ba = zeros(101,1);
vErrors_ba = pErrors_ba;

kk = 1;
for j = -5:0.1:5  
    PIM2 = PreintegrateMeasurement();
    PIM2.ba_ = j*ones(3,1);
    for i = 1:100
        
        PIM2 = PIM2.Preintegrate(acc_noise(i,:)',gyro_noise(i,:)',imuPara,step);
    end
    correctVMeas = PIM.v_ + PIM.Dv_ba_*PIM2.ba_;
    correctPMeas = PIM.p_ + PIM.Dp_ba_*PIM2.ba_;
    pErrors_ba(kk) = norm(correctPMeas-PIM2.p_);
    vErrors_ba(kk) = norm(correctVMeas-PIM2.v_);
    kk = kk + 1;
    % fprintf('linearization error: %f %f\n', norm(correctVMeas-PIM2.v_),norm(correctPMeas-PIM2.p_))
end
j = -5:0.1:5;
plot(j,pErrors_ba,j,vErrors_ba);
grid on
legend('position error w.r.t ba','velocity error w.r.t. ba')

pErrors_bg = zeros(101,1);
vErrors_bg = pErrors_bg;
phiErrors_bg = pErrors_bg;
kk = 1;
for j = -5:0.1:5  
    PIM2 = PreintegrateMeasurement();
    PIM2.bg_ = j*ones(3,1)*pi/180;
    for i = 1:100
        PIM2 = PIM2.Preintegrate(acc_noise(i,:)',gyro_noise(i,:)',imuPara,step);
    end
    correctRMeas = PIM.R_*SO3.exp(PIM.DR_bg_*PIM2.bg_);
    correctVMeas = PIM.v_ + PIM.Dv_bg_*PIM2.bg_;
    correctPMeas = PIM.p_ + PIM.Dp_bg_*PIM2.bg_;
    phiErrors_bg(kk) = norm(SO3.log(correctRMeas) - PIM2.phiv_);
    pErrors_bg(kk) = norm(correctPMeas - PIM2.p_);
    vErrors_bg(kk) = norm(correctVMeas - PIM2.v_);
    kk = kk + 1;
    % fprintf('linearization error: %f %f\n', norm(correctVMeas-PIM2.v_),norm(correctPMeas-PIM2.p_))
end
figure(2)
j = -5:0.1:5;
plot(j,pErrors_bg,j,vErrors_bg,j,phiErrors_bg);
grid on
legend('position error w.r.t bg','velocity error w.r.t. bg','oritation error w.r.t. bg')

% si =  NavState(zeros(3,1),r(1,:)',v(1,:)',zeros(3,1),zeros(3,1));
% sj = si.predict(PIM);
% % sj = si;
% cov = zeros(15,15);
% cov(1:9,1:9) = PIM.cov_;
% cov(10:15,10:15) = 0.1*eye(6);
% infor = eye(15)/cov;
% 
% oriMeas = SO3.log(AttitudeBase.a2cnb(atti(100,:))');
% for i = 1:10
%     [ err,J ] = IMUerror( si,sj,PIM,r(100,:)',oriMeas );
%     fprintf('err: %f\n',norm(err));
% %     if norm(err)<0.00001
% %         break;
% %     end
%     delta = -(J'*infor*J)\J'*infor*err;
% %     fprintf('delta: %f\n',norm(delta))
%     
%     sj = sj.update(delta);
%     si.ba_ = sj.ba_;
%     si.bg_ = sj.bg_;
% end



