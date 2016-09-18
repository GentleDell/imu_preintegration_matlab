classdef IMUPara
   properties
      accCov_;
      gyroCov_;
      bAccCov_;
      bGyroCov_;
      
   end
   
   methods
       function o = IMUPara(accCov,gyroCov,bAccCov,bGyroCov)
           o.accCov_ = accCov;
           o.gyroCov_ = gyroCov;
           o.bAccCov_ = bAccCov;
           o.bGyroCov_ = bGyroCov;
       end
       
   end
    
end