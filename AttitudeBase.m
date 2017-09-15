classdef AttitudeBase
    %������װ�˻����̬�仯��ϵ����
    
    methods(Static = true)
        %%
        % Problem?
        function cnb = a2cnb(atti)  % 3D rotation matrix
            sp = sin(atti(1));cp = cos(atti(1));    %pitch
            sr = sin(atti(2));cr = cos(atti(2));    %roll
            sy = sin(atti(3));cy = cos(atti(3));    %yaw
            cnb = zeros(3);
            cnb(1,1) = cr*cy - sr*sp*sy;
            cnb(1,2) = cr*sy + sr*sp*cy;
            cnb(1,3) = -sr*cp;
            cnb(2,1) = -cp*sy;
            cnb(2,2) = cp*cy;
            cnb(2,3) = sp;
            cnb(3,1) = sr*cy + cr*sp*sy;
            cnb(3,2) = sr*sy - cr*sp*cy;
            cnb(3,3) = cr*cp;
        end
        %%
        
        function a = cnb2atti(cnb )
            %a2cnb()�ķ�����
            gama = atan(-cnb(1,3)/cnb(3,3));
            theta = asin(cnb(2,3));
            ctheta = cos(theta);
            psi = acos(cnb(2,2)/ctheta);
            if -cnb(2,1)/ctheta < 0
                psi = -psi;
            end
            a = [theta gama psi];
            
        end
        %%
        
        function cnb = quat2cnb( quat )
            %���룺��̬��Ԫ��
            %���������ϵ����ϵ����̬�任����(navigation frame to body frame)
            cnb(1,1) = quat(1)*quat(1) + quat(2)*quat(2) - quat(3)*quat(3) - quat(4)*quat(4);
            cnb(1,2) = 2*(quat(2)*quat(3)+quat(1)*quat(4));
            cnb(1,3) = 2*(quat(2)*quat(4)-quat(1)*quat(3));
            cnb(2,1) = 2*(quat(2)*quat(3)-quat(1)*quat(4));
            cnb(2,2) = quat(1)*quat(1)-quat(2)*quat(2)+quat(3)*quat(3)-quat(4)*quat(4);
            cnb(2,3) = 2*(quat(3)*quat(4)+quat(1)*quat(2));
            cnb(3,1) = 2*(quat(2)*quat(4)+quat(1)*quat(3));
            cnb(3,2) = 2*(quat(3)*quat(4)-quat(1)*quat(2));
            cnb(3,3) = quat(1)*quat(1) - quat(2)*quat(2) - quat(3)*quat(3) + quat(4)*quat(4);
        end
        %%
        
        function A = w2Datti( atti )
            %���룺��̬��
            %�������̬����������ٶȵĹ�ϵ����
            % A = [-sin(atti(2)), 0, cos(atti(2));
            %     cos(atti(2))*cos(atti(1)),0,sin(atti(2))*cos(atti(1));
            %     sin(atti(1))*sin(atti(2)), cos(atti(1)), -sin(atti(1))*cos(atti(2))]/cos(atti(1));
            sp = sin(atti(1));cp = cos(atti(1));
            sr = sin(atti(2));cr = cos(atti(2));
            A = [cp*cr, 0, sr*cp;
                sp*sr, cp, -cr*sp;
                -sr, 0, cr]/cp;
        end
        %%
        % Problems?
        function T = Datti2w(atti)
            sp = sin(atti(1));cp = cos(atti(1));
            sr = sin(atti(2));cr = cos(atti(2));
            T = [cr, 0 -sr*cp;
                0, 1,sp;
                sr, 0, cp*cr];
            
        end
        %%
        
        function deltheta = w2omega( w )
            %���룺���ٶȵı仯��
            %�������Ԫ��仯������Ԫ��Ĺ�ϵ����
            delthetax = w(1);
            delthetay = w(2);
            delthetaz = w(3);
            deltheta(1,1) = 0;
            deltheta(1,2) = -delthetax;
            deltheta(1,3) = -delthetay;
            deltheta(1,4) = -delthetaz;
            deltheta(2,1) = delthetax;
            deltheta(2,2) = 0;
            deltheta(2,3) = delthetaz;
            deltheta(2,4) = -delthetay;
            deltheta(3,1) = delthetay;
            deltheta(3,2) = -delthetaz;
            deltheta(3,3) = 0;
            deltheta(3,4) = delthetax;
            deltheta(4,1) = delthetaz;
            deltheta(4,2) = delthetay;
            deltheta(4,3) = -delthetax;
            deltheta(4,4) = 0;
            
        end
        
    end
    
    
end
