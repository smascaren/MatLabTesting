function [x_out, t_out, U_out] = diffusion1d(kappa, x_rng, nx, t_rng, nt, u_init, u_bndry)
    % Unpack the space and time intervals
    a = x_rng(1);
    b = x_rng(2);
    t0 = t_rng(1);
    tfinal = t_rng(2);

    % Create space and time discretization
    x_out = linspace(a, b, nx);
    t_out = linspace(t0, tfinal, nt);
    % dx = x_out(2) - x_out(1);
    % dt = t_out(2) - t_out(1);

    dt = 0.0005;
    dx = 0.01;

    % Initialize the solution matrix
    U_out = zeros(nx, nt);

    % Set initial condition
    U_out(:,1) = u_init(x_out);

    % Calculate the stability criterion
    r = kappa * dt / dx^2;
    if r > 0.5
        error('The scheme is unstable with r = %.2f. Choose smaller dt or larger dx.', r);
    end

    % Time-stepping loop
    for n = 1:nt-1
        % Boundary conditions
        bndry_vals = u_bndry(t_out(n+1));
        U_out(1, n+1) = bndry_vals(1);
        U_out(end, n+1) = bndry_vals(2);

        % Update interior points using the modified explicit finite difference scheme
        for i = 2:nx-1
            U_out(i, n+1) = U_out(i, n) + (kappa * dt^2 / dx^3) * (U_out(i+1, n) - 2*U_out(i, n) + U_out(i-1, n));
        end
    end
end

% Initial condition function
function [u] = u2a_init(x)
    u = 0 * x + 0.9;
end

% Boundary condition function
function [u] = u2a_bndry(t)
    u = [1.7;
         4.1];
end

% % Parameters
% kappa = 0.1; % Example thermal diffusivity
% x_rng = [0, 1];
% nx = 6;
% t_rng = [1, 3];
% nt = 12;

% Call the diffusion function
% [x_out, t_out, U_out] = diffusion1d(kappa, x_rng, nx, t_rng, nt, @u2a_init, @u2a_bndry);
[x2, t2, U2] = diffusion1d( 0.1, [0,1], 100, [1, 3], 4000, @u2a_init, @u2a_bndry );


figure;
[X, T] = meshgrid(x_out, t_out);
mesh(T, X, U_out');
xlabel('Time');
ylabel('Space');
zlabel('Temperature');
title('Heat Diffusion in 1D');
































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Example usage
% kappa = 0.01;
% x_rng = [0, 10];
% nx = 100;
% t_rng = [0, 1];
% nt = 500;
% 
% [x_out, t_out, U_out] = diffusion1d(kappa, x_rng, nx, t_rng, nt, @u2a_init, @u2a_bndry);
% 
% % Plot the results (example plot at final time)
% figure;
% plot(x_out, U_out(:, end));
% title('Temperature distribution at final time');
% xlabel('Position (x)');
% ylabel('Temperature (u)');
% 
% 
% 
% 
% 
% function [x_out, t_out, U_out] = diffusion1d( kappa, x_rng, nx, t_rng, nt, u_init, u_bndry )
%     % Initialization
%     % ==============
%     %
%     % Set up the spatial and time grids, and initialize the solution matrix.
% 
%     % Define x_out and t_out
%     x_out = linspace(x_rng(1), x_rng(2), nx);
%     t_out = linspace(t_rng(1), t_rng(2), nt);
% 
%     % Initialize U_out
%     U_out = zeros(nx, nt);
% 
%     % Set initial condition
%     U_out(:, 1) = u_init(x_out)';
% 
%     % Calculate time step size (dt) and spatial step size (dx)
%     dt = (t_rng(2) - t_rng(1)) / (nt - 1);
%     dx = (x_rng(2) - x_rng(1)) / (nx - 1);
% 
%     % Calculate the stability parameter
%     r = kappa * dt / dx^2;
% 
%     if r > 0.5
%         warning('The scheme might be unstable. Consider reducing dt or increasing dx.');
%     end
% 
%     % Solving
%     % =======
%     %
%     % Use a finite difference method to solve the heat equation over time.
% 
%     % Time-stepping loop using explicit finite difference scheme
%     for j = 2:nt
%         % Apply boundary conditions at current time step
%         [u_bndry_left, u_bndry_right] = u_bndry(t_out(j));
%         U_out(1, j) = u_bndry_left;
%         U_out(end, j) = u_bndry_right;
% 
%         % Update the interior points
%         for i = 2:nx-1
%             U_out(i, j) = U_out(i, j-1) + r * (U_out(i-1, j-1) - 2*U_out(i, j-1) + U_out(i+1, j-1));
%         end
%     end
% end
% function [u] = u2a_bndry(t)
%     % Boundary conditions function
%     % Returns the boundary values at the ends of the spatial interval
%     u = [0*t + 1.7;  % Left boundary
%          0*t + 4.1]; % Right boundary
% end
% 
% function [u] = u2a_init(x)
%     % Initial condition function
%     % Returns the initial temperature distribution along the spatial interval
%     u = 0*x + 0.9;  % Initial temperature at each point in space
% end
