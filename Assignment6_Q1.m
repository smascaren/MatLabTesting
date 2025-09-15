% Sharlene M.
% Q1

clc
clear

% Given values
x = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50];
y = [17, 24, 31, 33, 37, 37, 40, 40, 42, 41];

% (a)
matr = ones(size(x')); % creating a matrix
trans = x'; %transpose
A = [trans matr];
b = y';
straight_line = A\b; % Using the formula

% Ploting for a straight line fit
figure(1)
subplot(2,1,1);
plot(x,straight_line(1)*x + straight_line(2))
hold on
title("Straight line");
plot(x,y, "o")
hold off

% (b)
f = @(fit,x) fit(1)*x.^fit(2); % Formula for the power equation
fit = lsqcurvefit(f,[1,1],x,y); % Using a built-in function

% Ploting for a power equation line
figure(2)
plot(x,f(fit,x))
hold on
title("Power Equation");
plot(x,y,"o")
hold off

% % (c)
f1 = @(fit1,x) fit1(1)*(x./(fit1(2)+x)); % Formula for the power equation
fit1 = lsqcurvefit(f1,[50 0.1],x,y); % Using a built-in function

% Ploting for a saturation-growthrate equation
figure(3)
x_val = linspace(0,50,100);
plot(x_val,f1(fit1,x_val));
hold on;
plot(x,y,'o');
title("Saturation-Growthrate Equation")
hold off

% (d)
fit2 = polyfit(x,y,2);
range = (0:55);
y_val = polyval(fit2,range); 

% Ploting for a parabola
figure(4)
plot(range, y_val);
hold on;
plot(x,y,'o');
title("Parabola")
hold off



