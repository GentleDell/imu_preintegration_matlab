clc;clear;
gt = KITTIGT.LoadPose('00.txt');
interN = 5;
N = size(gt,3);
interPose = IMURawData.SplineInterplotePose(gt,5);
step = 0.1;
[acc, gyro] = IMURawData.DifferentialIMUKinectics(interPose,step/interN);