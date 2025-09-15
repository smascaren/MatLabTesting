% Sharlene M.
% 2.

clc
clear

% All equations
matrix_A = [9 3 1; 2 5 -1; -6 0 8];
matrix_B = [4 2 -2; 1 5 -1; 1 1 6];
matrix_C = [-3 4 5; -2 2 -3; 0 2 -1];

% Get the diagonals from the matrix
x = abs(diag(matrix_A));
y = abs(diag(matrix_B));
z = abs(diag(matrix_C));


% For matrix A
for j = 1:3
    k = 0;
    iter = 0;
    for i = 1:3
     
        if i ~= j % If the matrix is diagonally dominant.
            k = k+abs(matrix_A(j,i)); 
        else
            iter = iter +1;

        end

    end
    
end
if x(j)>k % If the matrix is not diagonally dominant, thus not converging.
       disp("Set 1 system of equations converge.")
else
       disp("Set 1 system of equations do not converge.")
end


disp("********************************************")

% For matrix B
for j = 1:3
    k = 0;
    iter = 0;
    for i = 1:3
        
        if i ~= j % If the matrix is diagonally dominant.
            k = k+abs(matrix_B(j,i));
        else
            iter = iter +1;

        end

    end
    
end
if y(j)>k % If the matrix is not diagonally dominant, thus not converging.
       disp("Set 2 system of equations converge.")
else
       disp("Set 2 system of equations do not converge.")
end


disp("********************************************")

% For matrix C
x = 0;
y = 0;
z = 0;
max_iter =20;
iter =0;
while iter < max_iter % Subbing values into the equations to check if it converges later.
    x = x+1;
    y = x+1;
    z = z+1;
    f1 = (-3)*x+4*x+5*x-6;
    f2 = (-2)*y+2*y+(-3)*y+3;
    f3 = 2*z-z-1;

    if f1 == f2 && f2 == f3
        disp(iter)
        break
    else
        iter = iter +1;
    end
end
disp("Set 3 system of equations do not converge.") % if the system of equations does not converge.

