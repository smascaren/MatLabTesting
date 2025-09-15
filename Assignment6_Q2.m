% Sharlene M
% Q2

clc
clear

% Given values
x = [1.6, 2, 2.5, 3.2, 4, 4.5];
fx = [2, 8, 14, 15, 8, 2];

%first divided differences
fx2x1=(fx(2)-fx(1))/(x(2)-x(1));%this is b2
fx3x2=(fx(3)-fx(2))/(x(3)-x(2));
fx4x3=(fx(4)-fx(3))/(x(4)-x(3));

%second divided differences
fx3x2x1=(fx3x2-fx2x1)/(x(3)-x(1));%this is b3
fx4x3x2=(fx4x3-fx3x2)/(x(4)-x(2));

%third divided difference 
fx4x3x2x1=(fx4x3x2-fx3x2x1)/(x(4)-x(1));%this is b4


%newton interpol
b1=fx(1);
b2=fx2x1;
b3=fx3x2x1;
b4=fx4x3x2x1;

% Getting the value for each order and displaying it with the value of 2.8.
f_1=@(X) b1+b2.*(X-x(1));
f_2=@(X) b1+b2.*(X-x(1))+b3.*(X-x(1)).*(X-x(2));
f_3=@(X) b1+b2.*(X-x(1))+b3.*(X-x(1)).*(X-x(2))+b4.*(X-x(1)).*(X-x(2)).*(X-x(3));

disp("First order: " +f_1(2.8))
disp("Second order: " +f_2(2.8))
disp("Third order: " +f_3(2.8))


%******************************************************************
% Lagrange interpolation

%Formula with the fourth points to get the third order.
f2_lagrange=@(X) (((X-x(2)).*(X-x(3)).*(X-x(4)))/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))))*fx(1)+...
                 (((X-x(1)).*(X-x(3)).*(X-x(4)))/((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))))*fx(2)+...
                 (((X-x(1)).*(X-x(2)).*(X-x(4)))/((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))))*fx(3)+...
                 (((X-x(1)).*(X-x(2)).*(X-x(3)))/((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))))*fx(4);

% Displaying the results.
disp(" ")
disp("Lagrange interpolation 3rd order:" + f2_lagrange(2.8))

% Stating the conclusion drawn from these results.
disp("Yes the values for the Newton and Lagrange interpolation to the third order is the same.")


