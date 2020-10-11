%The function calculates the first order derivative in a specific direction 
%'dx' at a specific point 'J' given a matrix 'y'  
function dydx=diffmy(J,y,dx)
n=length(y);
   if J==1
           dydx=(-y(3)+4*y(2)-3*y(1))/(2*dx);
           %dydx=(y(J+1)-y(n-1))/(2*dx);
   else
       if J==n
           dydx=(3*y(J)-4*y(J-1)+y(J-2))/(2*dx);
           %dydx=(y(2)-y(J-1))/(2*dx);
       else
           dydx=(y(J+1)-y(J-1))/(2*dx);
       end
   end
end
