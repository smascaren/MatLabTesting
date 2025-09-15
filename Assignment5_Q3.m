% Sharlene M.
% 3.

clc
clear

% (b)

% System 1

x = -10:0.1:10;
xlim([-10 10])
f1 = -x.^2 + 10;
f2 = -x + 15;

% Plot the equations
figure (1)
plot(x,f1) % plot the parabola
hold on
plot(x,f2) % plot the line
hold off
title('System 1: Parabola and a Line');

disp("There is no solution for System 1 has both functions do not intersect at any points.")

% System 2

f3 = @(x,y) y^2 + x^2 - 26;
f4 = @(x,y) 25*y^2 + 3*x^2 - 100;

% Plot the equations
figure(2)
ezplot(f3,[-6,6,-6,6]);  % plot ellipse
hold on;
ezplot(f4,[-6,6,-6,6]);  % plot circle
title('System 2: Ellipse and Circle');

x1 = [1,1];
x2 = [-1,-1]; 
x3 = [-1,1];
x4 = [1,-1];

disp("These system of equations have 4 solutions:")

% Solve the system of equations
sol1 = fsolve(@(x) [f3(x(1),x(2)); f4(x(1),x(2))], x1)
sol2 = fsolve(@(x) [f3(x(1),x(2)); f4(x(1),x(2))], x2)
sol3 = fsolve(@(x) [f3(x(1),x(2)); f4(x(1),x(2))], x3)
sol4 = fsolve(@(x) [f3(x(1),x(2)); f4(x(1),x(2))], x4)

