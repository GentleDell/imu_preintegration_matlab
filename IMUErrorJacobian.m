classdef IMUErrorJacobian
    
    methods(Static = true)
        function [ err,J ] = GPSAndOri( pi,pj,imuMeas,gpsMeas,oriMeas )
            err = zeros(15,1);
            J = zeros(15,15);
            g = [0,0,-9.8]';
            delta_bg = pi.bg_ - imuMeas.bg_;
            delta_ba = pi.ba_ - imuMeas.ba_;
            %     err(10:15) = [delta_ba;delta_bg];
            
            t = imuMeas.t_;
            t22 = t*t/2;
            
            phiDelta_bg = imuMeas.DR_bg_*delta_bg;
            correctRMeas = imuMeas.R_*SO3.exp(phiDelta_bg);
            correctVMeas = imuMeas.v_ + imuMeas.Dv_ba_*delta_ba + imuMeas.Dv_bg_*delta_bg;
            correctPMeas = imuMeas.p_ + imuMeas.Dp_ba_*delta_ba + imuMeas.Dp_bg_*delta_bg;
            
            err(1:3) = SO3.log(correctRMeas'*pi.R_'*pj.R_);
            err(4:6) = pi.R_'*(pj.v_ - pi.v_ - g*t) - correctVMeas;
            err(7:9) = pi.R_'*(pj.p_ - pi.p_ - pi.v_*t - g*t22) - correctPMeas;
            err(10:12) = pj.p_ - gpsMeas;
            err(13:15) = pj.phiv_ - oriMeas;
            
            Jr_Err_errPhiInv = SO3.Dlog(err(1:3));
            RErrPhi = SO3.exp(err(1:3));
            Jr_phiDelta_bg = SO3.Dexp(phiDelta_bg);
            J(1:3,1:3) = Jr_Err_errPhiInv;
            J(1:3,13:15) = -Jr_Err_errPhiInv*RErrPhi'*Jr_phiDelta_bg*imuMeas.DR_bg_;
            
            J(4:6,4:6) = pi.R_';
            J(4:6,10:12) = -imuMeas.Dv_ba_;
            J(4:6,13:15) = -imuMeas.Dv_bg_;
            
            %     J(7:9,7:9) = pi.R_'*pj.R_;
            J(7:9,7:9) = pi.R_';
            J(7:9,10:12) = -imuMeas.Dp_ba_;
            J(7:9,13:15) = -imuMeas.Dp_bg_;
            
            J(10:12,7:9) = eye(3);
            J(13:15,1:3) = eye(3);
        end
        
        function [ err,J ] = GPS( pi,pj,imuMeas,gpsMeas )
            err = zeros(12,1);
            J = zeros(12,15);
            g = [0,0,-9.8]';
            delta_bg = pi.bg_ - imuMeas.bg_;
            delta_ba = pi.ba_ - imuMeas.ba_;
            %     err(10:15) = [delta_ba;delta_bg];
            
            t = imuMeas.t_;
            t22 = t*t/2;
            
            % incorporating Bias update
            phiDelta_bg = imuMeas.DR_bg_*delta_bg;
            correctRMeas = imuMeas.R_*SO3.exp(phiDelta_bg);
            correctVMeas = imuMeas.v_ + imuMeas.Dv_ba_*delta_ba + imuMeas.Dv_bg_*delta_bg;
            correctPMeas = imuMeas.p_ + imuMeas.Dp_ba_*delta_ba + imuMeas.Dp_bg_*delta_bg;
            
            % residual error besed on formula (37) in the paper
            % [R, v, p, GNSS]'
            err(1:3) = SO3.log(correctRMeas'*pi.R_'*pj.R_);
            err(4:6) = pi.R_'*(pj.v_ - pi.v_ - g*t) - correctVMeas;
            err(7:9) = pi.R_'*(pj.p_ - pi.p_ - pi.v_*t - g*t22) - correctPMeas;     % g need 1/2 ?
            err(10:12) = pj.p_ - gpsMeas;
            
            % obtain Jaccobian of residual error
            % state order is [fai, pose, velovcity], derivative variances
            % is [Fia_i, pi, vi, Fia_j, pj, vj, ba, bg]
            Jr_Err_errPhiInv = SO3.Dlog(err(1:3));
            RErrPhi = SO3.exp(err(1:3));
            Jr_phiDelta_bg = SO3.Dexp(phiDelta_bg);
            J(1:3,1:3) = Jr_Err_errPhiInv;          % d resierror_deltaRij, d Fiaj
            J(1:3,13:15) = -Jr_Err_errPhiInv*RErrPhi'*Jr_phiDelta_bg*imuMeas.DR_bg_;    % d resierror_deltaRij, d bg
            
            J(4:6,4:6) = pi.R_';    % d resierror_deltavij, d vj
            J(4:6,10:12) = -imuMeas.Dv_ba_;     % d resierror_deltavij, d ba
            J(4:6,13:15) = -imuMeas.Dv_bg_;     % d resierror_deltavij, d bg
            
            %     J(7:9,7:9) = pi.R_'*pj.R_;    % ??? why not use this
            J(7:9,7:9) = pi.R_';
            J(7:9,10:12) = -imuMeas.Dp_ba_;     % d resierror_deltavij, d ba
            J(7:9,13:15) = -imuMeas.Dp_bg_;     % d resierror_deltavij, d bg
            
            J(10:12,7:9) = eye(3);              % ?? minus
        end
        
    end
    
end



