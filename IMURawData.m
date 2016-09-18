classdef IMURawData
    
    methods
        
        
    end
    
    methods(Static = true)
        
        function [acc, gyro] = DifferentialIMUKinectics(groundtruth,step)
            N = size(groundtruth,3);
            acc = zeros(N - 2,3);
            gyro = acc;
            g = [0 0 -9.8]';
            for i = 1:N-2
                deltaR = groundtruth(1:3,1:3,i)'*groundtruth(1:3,1:3,i+1);
                gyro(i,:) = SO3.log(deltaR)'/step;
                a = (groundtruth(:,4,i+2) - 2*groundtruth(:,4,i+1) + groundtruth(:,4,i))/step/step;
                acc(i,:) = (a-g)'*groundtruth(1:3,1:3,i);
            end
        end
        
        function interploteT = SplineInterplotePose(groundtruth,interN)
            N = size(groundtruth,3);
            interNinv = 1/interN;
            interNinv2 = interNinv*interNinv;
            interNinv3 = interNinv2*interNinv;
            
            Bu = zeros(interN, 4);
            C = [6 0 0 0;
                5 3 -3 1;
                1 3 3 -2;
                0 0 0 1]/6;
            for j = 1:interN-1
                Bu(j,:) = [1 j*interNinv j^2*interNinv2 j^3*interNinv3]*C';
            end
            interploteT = zeros(3,4,N-4+(N-5)*(interN-1));
            
            k = 1;
            for i = 1:N-4
                
                if i == 1
%                     T_i = [eye(3) zeros(3,1)];
                    T_i = groundtruth(:,:,i);
                    T_i_1 = groundtruth(:,:,i);
                    T_i_2 = groundtruth(:,:,i+1);
                    T_i_3 = groundtruth(:,:,i+2);
                    
                    omega_i = SE3.log(SE3.multiply(SE3.inv(T_i),T_i_1));
                    omega_i_1 = SE3.log(SE3.multiply(SE3.inv(T_i_1),T_i_2));
                    omega_i_2 = SE3.log(SE3.multiply(SE3.inv(T_i_2),T_i_3));
                else
                    T_i = groundtruth(:,:,i-1);
                    omega_i = omega_i_1;
                    omega_i_1 = omega_i_2;
                    omega_i_2 = SE3.log(SE3.multiply(SE3.inv(groundtruth(:,:,i+1)),groundtruth(:,:,i+2)));
                end
                interploteT(:,:,k) = groundtruth(:,:,i);
                k = k + 1;
                for j = 1:interN-1
                    T1 = SE3.exp(Bu(j,2)*omega_i);
                    T2 = SE3.exp(Bu(j,3)*omega_i_1);
                    T3 = SE3.exp(Bu(j,4)*omega_i_2);
                    interploteT(:,:,k) = SE3.multiply(T_i,SE3.multiply(T1,SE3.multiply(T2,T3)));
                    k = k + 1;
                end
                
            end
        end
    end
    
end