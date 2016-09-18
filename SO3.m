classdef SO3
    properties
        
    end
    
    methods(Static = true)
        function R = exp(phiv)
            phi = norm(phiv);
            phix = SO3.skew(phiv);
            phiInv = 1/phi;
            I3 = eye(3);
            if(phi<1e-3)
                R = I3 + phix;
            else
                R = I3 + sin(phi)*phiInv*phix+(1-cos(phi))*phiInv*phiInv*phix*phix;
            end
        end
        
        function phiv = log(R)
            tr = trace(R);
            if(abs(tr+1)<1e-10)
                if(abs(R(3,3)+1)>1e-10)
                    phiv = (pi/sqrt(2+2*R(3,3)))*[R(1,3),R(2,3),1+R(3,3)]';
                else if(abs(R(2,2)+1)>1e-10)
                        phiv = (pi/sqrt(2+2*R(2,2)))*[R(1,2),1+R(2,2),R(3,2)]';
                    else
                        phiv = (pi/sqrt(2+2*R(1,1)))*[1+R(1,1),R(2,1),R(3,1)]';
                    end
                end
            else
                tr_3 = tr - 3.;
                if(tr_3<-1e-7)
                    theta = acos((tr-1)/2);
                    magnitude = theta/(2*sin(theta));
                else
                    magnitude = 0.5 - tr_3*tr_3/12;
                end
                phiv = magnitude*[R(3,2)-R(2,3),R(1,3)-R(3,1),R(2,1)-R(1,2)]';
            end
        end
        
        function phix = skew(phi)
            phix = [0,-phi(3),phi(2);
                phi(3),0,-phi(1);
                -phi(2),phi(1),0];
        end
        
        function Jr = Dexp(phiv)
            phi = norm(phiv);
            I3 = eye(3);
            phix = SO3.skew(phiv);
            phiinv = 1/phi;
            if(phi<1e-3)
                Jr = I3 - 0.5*phix;
            else
                Jr = I3 - (1-cos(phi))*phiinv*phiinv*phix+(phi-sin(phi))*phiinv*phiinv*phiinv*phix*phix;
            end
        end
        
        function JrInv = Dlog(phiv)
            phi = norm(phiv);
            I3 = eye(3);
            phix = SO3.skew(phiv);
            phiinv = 1/phi;
            if(phi<1e-3)
                JrInv = I3 + 0.5*phix;
            else
                JrInv = I3 + 0.5*phix + (phiinv*phiinv-(1+cos(phi))*phiinv/(2*sin(phi)))*phix*phix;
            end
        end
        
        function Jl = DexpLeft(phiv)
            phi = norm(phiv);
            phiv = phiv/phi;
            I3 = eye(3);
            phix = SO3.skew(phiv);
            phiinv = 1/phi;
            if(phi<1e-3)
                Jl = I3;
            else
                a = sin(phi)/phi;
                Jl = a*I3+(1 - a)*(phiv*phiv') + (1 - cos(phi))*phiinv*phix;
            end
            
        end
        
        function JlInv = DlogLeft(phiv)
            
            phi = norm(phiv);
            phi2 = phi/2;
            phiv = phiv/phi;
            I3 = eye(3);
            phix = SO3.skew(phiv);
            if(phi2<1e-3)
                JlInv = I3;
            else
                JlInv = phi2*cot(phi2)*I3 + (1 - phi2*cot(phi2))*(phiv*phiv') - phix*phi2;
            end
            
        end
    end
    
end
