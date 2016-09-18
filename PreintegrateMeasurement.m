classdef PreintegrateMeasurement
    properties
        R_;
        p_;
        v_;
        phiv_;
        ba_;
        bg_;
        cov_;
        t_;
        DR_bg_;
        Dv_bg_;
        Dp_bg_;
        Dv_ba_;
        Dp_ba_;
    end
    
    methods
        function o = PreintegrateMeasurement()
            o.R_ = eye(3);
            o.p_ = zeros(3,1);
            o.v_ = zeros(3,1);
            o.phiv_ = zeros(3,1);
            o.ba_ = o.phiv_;
            o.bg_ = o.phiv_;
            o.cov_ = eye(9);
            o.t_ = 0;
            o.DR_bg_ = zeros(3);
            o.Dv_bg_ = o.DR_bg_;
            o.Dp_bg_ = o.DR_bg_;
            o.Dv_ba_ = o.DR_bg_;
            o.Dp_ba_ = o.DR_bg_;
        end
        
        function o = Preintegrate(o,acc,gyro,imuPara,t)
            o.t_ = o.t_ + t;
            I3 = eye(3);
            z3 = zeros(3);
            t22 = t*t/2;
            %correct gyro measurement with bias
            gyroCorrect = gyro - o.bg_;
            %calculate the orientation change between consecutive time slots
            Rk_k_1 = SO3.exp(gyroCorrect*t);
            Rk_k_1T = Rk_k_1';
            %calculate Jacobian of Jr
            Jrk_k_1 = SO3.Dexp(gyroCorrect*t);
            %correct acc measurement with bias
            accCorrect = acc - o.ba_;
            accCorrectX = SO3.skew(accCorrect);
            
            %update state covariance with new IMU measurements
            Ak = [Rk_k_1T,z3,z3;
                -o.R_*accCorrectX*t22,I3,I3*t;
                -o.R_*accCorrectX*t,z3,I3];
            Bk = [z3;
                o.R_*t22;
                o.R_*t];
            Ck = [Jrk_k_1;
                z3;
                z3];
            o.cov_ = Ak*o.cov_*Ak' + Bk*imuPara.accCov_*Bk' + Ck*imuPara.gyroCov_*Ck';
            
            DR_bg = o.DR_bg_;
            Dv_bg = o.Dv_bg_;
            Dp_bg = o.Dp_bg_;
            Dv_ba = o.Dv_ba_;
            Dp_ba = o.Dp_ba_;
            
            o.DR_bg_ = Rk_k_1T*DR_bg - Jrk_k_1*t;
            o.Dv_bg_ = Dv_bg - o.R_*accCorrectX*DR_bg*t;
            o.Dv_ba_ = Dv_ba - o.R_*t;
            o.Dp_bg_ = Dp_bg + Dv_bg*t - o.R_*accCorrectX*DR_bg*t22;
            o.Dp_ba_ = Dp_ba + Dv_ba*t - o.R_*t22;
            
            %update navigation states with IMU measurements
            R = o.R_;
            o.R_ = o.R_*Rk_k_1;
            o.phiv_ = SO3.log(o.R_);
            o.p_ = o.p_ + o.v_*t + R*accCorrect*t22;
            o.v_ = o.v_ + R*accCorrect*t;
        end
        
    end
    
end