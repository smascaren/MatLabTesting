% Sharlene Mascarenhas (21011314)
% 1.

clc
clear

maxit = 50; % maximum number of iterations
tol = 1e-6;
eps = 1;  % initialize epsilon for error, any value OK as long as greater than the tolerance

%(a)
f = @(x) exp(-0.5*x)*(4-x)-2;   % define function
fderiv = @(x) 0.5*exp(-0.5*x)*x-(3*exp(-0.5*x));   % define derivative of f

x=zeros(1,maxit);

%(b) initial guesses - Observation about new guesses: The bigger the guess the less iterations it
%has to take to find the roots.


x(1) = 2; % initial guess
% x(1) = 4; % initial guess
% x(1) = 6; % initial guess

i = 0;
while eps > tol && i <= maxit  % method for checking and getting the roots.
    i = i+1;
    x(i+1) = x(i) - f(x(i))/fderiv(x(i));
    eps = abs((x(i+1)-x(i))/x(i+1));
end
if eps > tol
    disp('Maximum number or iterations reached, cannot find any root')
else
    fprintf('For x0 = %1.6f, one root is %1.6f. \rThe number of iterations is %1.0f.\r',x(1),x(i),i)
end

dis = 1.7^2-4*(-0.9)*2.5;
rootcalc1 = (-1.7+sqrt(dis))/2/(-0.9);
rootcalc2 = (-1.7-sqrt(dis))/2/(-0.9);

fprintf('\rThe first calculated root is %1.6f \r',rootcalc1)
fprintf('\rThe second calculated root is %1.6f \r',rootcalc2)
