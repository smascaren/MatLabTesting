% Sharlene Mascarenhas (21011314)
% 2.3

clc
clear

% (d)
built_in = fsolve(@(x) myfun_conc(x),[1;1]);
disp("Using fsolve, the value is: " + built_in)

newton_conc(@(x) myfun_conc(x),@(x) myjac_conc((@(x) myfun_conc(x)),x),[1;1])

% (a)
function F = myfun_conc(x)
    %the input is a vector now 
	f1 = (4*10^-4)*(50-2*x(1))^2*(20-x(1))-(5+x(1)+x(2));
    f2 = (3.7*10^-2)*(50-2*x(1))*(10-x(2))-(5+x(1)+x(2));
	F = [f1;f2];

end

% (b)
function J1 = myjac_conc(F,x)
    new_x1 = x +[(1*10^-6)*x(1);0];
    new_x2 = x + [0;(1*10^-6)*x(2)];
    J1 = zeros(2);
    J1(:,1) = (F(new_x1)-F(x))/((1*10^-6)*x(1));

    J1(:,2) = (F(new_x2)-F(x))/((1*10^-6)*x(2));
  

end


% (c)
function [sol,iter] = newton_conc(f,J1,x0)
    % Solve a system of nonlinear equations f(x) = 0, 
    %Jac is the Jacobian 
    % x0 is a vector of the initial guesses
    maxiter = 10; % maximum number of iterations before exiting the program
    tol = 0.001; % tolerance for convergence
    eps = 1; % initial epsilon, any value OK as long as greater than the tolerance
    xold = x0; % initial guess
    i = 0; % initialize the number of iterations
 
    while eps > tol && i <= maxiter
        i = i+1; % increase the number of iterations
        Ja = J1(xold); % calculate the Jacobian at {xi} (i.e. xold)
        xnew = xold - inv(Ja)*f(xold); % calculate {xi+1} (i.e. xnew)
        eps = max(abs((xnew-xold)./xnew)); % determine the updated epsilon
        xold = xnew;
    end
    if eps > tol
        disp('Maxium number of iterations reached, cannot find solution')
        disp('Change the initial guesses')
        sol = [];
        iter = maxiter;
    else
        sol = xnew;
        iter = i;
        disp("The number of iteration completed: " + iter)
    end
end