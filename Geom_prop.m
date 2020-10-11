classdef Geom_prop < handle
    %Class of the geometry properties of a bearing
    properties %(Access=protected)%(GetAccess=public,SetAccess=public) 
        X1
        X2
        e_r
        r
        R
        h0
        L
        m
        x1
        x3
        s
    end
           
    methods %(Access=public)
        
        %Construct an instance of this class 
        function geom=Geom_prop(X1,X2,e_r,r,R,h0,L,m,x1,x3,s)
            if nargin>0
                geom.X1=X1;
                geom.X2=X2;
                geom.e_r=e_r;
                geom.r=r;
                geom.R=R;
                geom.h0=h0;
                geom.L=L;
                geom.m=m;
                geom.x1=x1;
                geom.x3=x3;
                geom.s=s;
            end
        end
                                  
        %Gap function
        function gap= fun_gap(geom)
            [h,hpow3,dhpow3dx1]=geomFunc(geom);
            %Results 
            gap.h=h;
            gap.hpow3=hpow3;
            gap.dhpow3dx1=dhpow3dx1;
        end
    end
end
                