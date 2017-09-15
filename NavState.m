classdef NavState
    properties
       R_;
       p_;
       v_;
       phiv_;
       bg_;
       ba_;
    end
    
    methods
        function o = NavState(phiv,p,v,bg,ba)
            o.R_ = SO3.exp(phiv);
            o.phiv_ = phiv;
            o.p_ = p;
            o.v_ = v;
            o.bg_ = bg;
            o.ba_ = ba;
        end
        
              
        function o = update(o,delta)
           R = o.R_;
           o.R_ = o.R_*SO3.exp(delta(1:3));
           o.phiv_ = SO3.log(o.R_);
%            o.p_ = o.p_ + R*delta(7:9);
           o.p_ = o.p_ + delta(7:9); 
           o.v_ = o.v_ + delta(4:6);
           o.ba_ = o.ba_ + delta(10:12);
           o.bg_ = o.bg_ + delta(13:15);
        end
        
        function pj = predict(o,PIM)
           g = [0;0;-9.8];
           t22 = PIM.t_*PIM.t_/2;
           % formula (31) of the paper
           R =  o.R_*PIM.R_;
           v = o.R_*PIM.v_+o.v_ + g*PIM.t_; 
           p = o.p_ + o.v_*PIM.t_ + g*t22 + o.R_*PIM.p_;
           pj = NavState(SO3.log(R),p,v,o.bg_,o.ba_);
           
        end
        
    end
    
end