% Define the grid size
n_x = 201;
n_y = 101;

% Initialize the disturbance
U6c_init = zeros(n_x, n_y);
dU6c_init = zeros(n_x, n_y);

% Set up the initial disturbance
for j = 1:n_x
    for k = 1:n_y
        x = 2 * (j - 1) / (n_x - 1);
        y = (k - 1) / (n_y - 1);
        
        % Initial disturbance centered at (0.5, 0.5)
        U6c_init(j, k) = exp(-1000 * ((x - 0.5)^2 + (y - 0.5)^2));
    end
end

% Run the simulation
[t6c, U6c] = wave2d(1, 1, U6c_init, dU6c_init, @U6c_bndry, [0, 350], 710);

% Function to set boundary conditions
function [U] = U6c_bndry(t, n_x, n_y)
    U = -Inf * ones(n_x, n_y);
 
    for j = 1:n_x
        for k = 1:n_y
            x = 2 * (j - 1) / (n_x - 1);
            y = (k - 1) / (n_y - 1);
        
            % Boundary conditions outside a circle of radius 0.5
            if sqrt((x - 1)^2 + (y - 0.5)^2) >= 0.5
                U(j, k) = NaN;
            end
        end
    end
end

% Function to compute the Laplacian
function L = laplacian(U, dx)
    L = (circshift(U, [1, 0]) + circshift(U, [-1, 0]) + ...
         circshift(U, [0, 1]) + circshift(U, [0, -1]) - ...
         4 * U) / dx^2;
end

% Wave propagation function
function [t, U_out] = wave2d(c, dx, U_init, dU_init, bndry_func, tspan, n_t)
    dt = (tspan(2) - tspan(1)) / (n_t - 1);
    t = linspace(tspan(1), tspan(2), n_t);
    
    [n_x, n_y] = size(U_init);
    U = U_init;
    dU = dU_init;
    U_next = U;
    
    U_out = zeros(n_x, n_y, n_t);
    U_out(:, :, 1) = U;
    
    for k = 2:n_t
        % Update boundary conditions
        U_bndry = bndry_func(t(k), n_x, n_y);
        U(U_bndry == NaN) = NaN;
        
        % Update wave equation
        L = laplacian(U, dx);
        U_next = 2 * U - dU + c^2 * dt^2 * L;
        
        % Enforce boundary conditions again
        U_next(U_bndry == NaN) = NaN;
        
        % Update variables for next time step
        dU = U;
        U = U_next;
        
        U_out(:, :, k) = U;
    end
end

% Find the time slice with the largest negative value at the point (1.5, 0.5)
target_x_index = round(1.5 * (n_x - 1) / 2) + 1;
target_y_index = round(0.5 * (n_y - 1)) + 1;

[z_min, t_min] = min(squeeze(U6c(target_x_index, target_y_index, :)));
disp(['The time slice with the largest negative value at (1.5, 0.5) is at t = ', num2str(t6c(t_min))]);

% Plot the results
figure;
surf(U6c(:,:,t_min));
title(['Wave Signal at Time Slice t = ', num2str(t6c(t_min))]);
xlabel('X Index');
ylabel('Y Index');
zlabel('Wave Amplitude');
grid on;
view(3);

% Animate the frames of the wave propagation
function frames = animate(U_out)
    [n_x, n_y, n_t] = size(U_out);
    frames = struct('cdata',[],'colormap',[]);
    
    figure;
    for k = 1:n_t
        surf(U_out(:,:,k));
        title(['Smascare']);
        xlabel('X ');
        ylabel('Y ');
        zlabel('Wave Amplitude');
        grid on;
        view(3);
        
        drawnow;
        frames(k) = getframe(gcf);
    end
end

% Display the animation
frames6c = animate(U6c);
