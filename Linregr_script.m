% Sharlene Mascarenhas (21011314)
% 2.1

clc
clear

% x and y values
x=[10 20 30 40 50 60 70 80];
y=[25 70 380 550 610 1220 830 1450];
[a, r2, syx] = linregr2(x,y)

% plots 

d = (a(1)+(a(2))*(x));
rs= (y-d);

subplot(2,1,1);
plot(x,d)
title("Subplot 1");
hold on
plot(x,y)
subplot(2,1,2);
plot(x,rs)
title("Subplot 2");
hold off

% Function for the method
function [a, r2, syx] = linregr2(x,y)

n=length(y); %Gets the length of the y values.

% Both A0 and a1 formulas
a_one = (n.*sum(x.*y)-(sum(x)).*(sum(y)))/(n.*(sum(x.^2))-(sum(x)).^2);
a_zero = (sum(y)/n) - a_one*(sum(x)/n)
h = y - (a_one*x)

% R formula
r = (n*sum(x.*y)-(sum(x)).*(sum(y)))/((sqrt(n.*sum(x.^2)-(sum(x)).^2)).*(sqrt(n.*sum(y.^2)-((sum(y)).^2))));
r2 = r.^2;

q = (y-a_zero-(a_one)*(x)).^2;
sr = sum(q);
syx = sqrt(sr/(n-2));

% Putting the a0 and a1 into one variable
a = [a_zero,a_one]
end


