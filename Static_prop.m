classdef Static_prop < Kinem_prop & handle
    %Class of the dynamic properties of a bearing
    
    properties(GetAccess=public,SetAccess=public)
        p0
        mu
        vel
    end
    
    methods
        %Construct an instance of this class 
        function static=...
                Static_prop(p0,mu,vel,n_rps,V1,V2,gap,X1,X2,...
                            e_r,r,R,h0,L,m,x1,x3,s)
                static=static@Kinem_prop(n_rps,V1,V2,gap,X1,X2,...
                            e_r,r,R,h0,L,m,x1,x3,s);
                if nargin>0
                    static.p0=p0;
                    static.mu=mu;
                    static.vel=vel;
                end
        end
 
        %Pressure field, total force and torque
        function pf=fun_pf(static)
            [p,F1,F2,M,Lambda,Gamma]=pressureField(static);
                
            %Results
            pf.p=p;
            pf.F1=F1;
            pf.F2=F2;
            pf.M=M;
            pf.Lambda=Lambda;
            pf.Gamma=Gamma;
            end
    end
end

        

            



