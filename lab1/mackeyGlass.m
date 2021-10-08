function [x] = mackeyGlass(beta,gamma,n,tau,tend)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x(1) = 1.5;

for t = 2:tend
   xold = x(t-1);
   
   if t-tau < 1
       x_tau = 0;
   else
       x_tau = x(t-tau);
   end
   x(t) = xold + (beta*x_tau)/(1+x_tau^n)-gamma*xold;
   
end


end

