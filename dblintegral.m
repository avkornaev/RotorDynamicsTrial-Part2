%The function calculates the double integral using Simpson's formula
%The integration limits are uniform from 0 to 1. That's why the result
%should be multiplied by the values of difference between the limits.
function [integ]=dblintegral(f)
%int=0;
[nx, ny]=size(f);
x1=linspace(0,1,nx);
y1=linspace(0,1,ny);
dx1=x1(2)-x1(1);
dy1=y1(2)-y1(1);
integ=0;
for i=2:2:nx-1
    for j=2:2:ny-1
         integ=integ+...
           (f(i-1,j-1)+f(i+1,j-1)+f(i-1,j+1)+f(i+1,j+1)+4*(f(i,j+1)+...
            f(i,j-1)+f(i-1,j)+f(i+1,j))+16*f(i,j))*dx1*dy1/9;
    end
%     if mod(nx,2)==0
%         integ=integ+
end

end
