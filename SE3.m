classdef SE3
   
   methods(Static = true)
       function T = exp(phiv)
           T = zeros(3,4);
           T(1:3,1:3) = SO3.exp(phiv(1:3));
           Jl = SO3.DexpLeft(phiv(1:3));
           T(:,4) = Jl*phiv(4:6);
       end
       
       function phiv = log(T)
           phiv = zeros(6,1);
           phiv(1:3) = SO3.log(T(1:3,1:3));
           JlInv = SO3.DlogLeft(phiv(1:3));
           phiv(4:6) = JlInv*T(:,4);
       end
       
       function poseInv = inv(pose)
           poseInv = [pose(1:3,1:3)' -pose(1:3,1:3)'*pose(:,4)];
       end
       
       function pose3 = multiply(pose2,pose1)
           pose3 = [pose2(1:3,1:3)*pose1(1:3,1:3) pose2(1:3,1:3)*pose1(:,4)+pose2(:,4)];
       end
       
   end
    
    
end