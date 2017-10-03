clc;clear;

groundtruth = KITTIGT.LoadPose('00.txt');
interN = 5; % but donot relate to time
N = size(groundtruth,3);
interPose = IMURawData.SplineInterplotePose(groundtruth,interN);     % 插值用于生成 IMU 数据
step = 0.1;
[acc, gyro] = IMURawData.DifferentialIMUKinectics(interPose,step/interN);

%% combine imu preintegration