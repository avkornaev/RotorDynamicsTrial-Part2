classdef Kinem_prop < Geom_prop & handle
    %Class of the kinematic properties of a bearing
    
    properties(GetAccess=public,SetAccess=public)
        n_rps
        V1
        V2
        gap
    end

    methods
         %Construct an instance of this class 
        function kinem=...
             Kinem_prop(n_rps,V1,V2,gap,X1,X2,e_r,r,R,h0,L,m,x1,x3,s)
                kinem=kinem@Geom_prop(X1,X2,e_r,r,R,h0,L,m,x1,x3,s);
                if nargin>0
                    kinem.n_rps=n_rps;
                    kinem.V1=V1;
                    kinem.V2=V2;
                    kinem.gap=gap;
                end
        end
        %Velocity functions  
        function vel=fun_vel(kinem)%Fill Mesh
            [u1,u2,du1hdx1] = velFunc(kinem);
            %Results
            vel.u1=u1;
            vel.u2=u2;
            vel.du1hdx1=du1hdx1;
        end
    end
end
                

