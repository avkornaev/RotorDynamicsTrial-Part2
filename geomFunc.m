function [h,hpow3,dhpow3dx1]=geomFunc(geom)
% Gap function and its derivatives

%Initialisation
X1=geom.X1;
X2=geom.X2;
x1=geom.x1;%x1 excluding x1=2*pi,as well as p(2*pi) is a known value
r=geom.r;
m=geom.m;
h0=geom.h0;
s=geom.s;

%Calculation
h=h0-X1*sin(x1/r)-X2*cos(x1/r);% gap function
hpow3=h.^3;% gap function powered by 3

dhpow3dx1=zeros(1,m(1));
for j=1:length(x1)
    dhpow3dx1(j)=diffmy(j,hpow3,s(1));% derivative d(h^3)/dx1
end

end