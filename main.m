% These control the number of iterations in both space and time. Note that
% this is not the number of steps, which would be either of these values
% minus 1. The output matrices will be exactly these dimensions.
t_iter_heat = 5000;
x_iter_heat = 75;

% The t-values used for analysis are contained in t_vals, and the x-values
% used for analysis are contained in x_vals. X_analytical contains the
% computed values for the analytical solution. The first row contains the
% initial condition, and each subsequent row contains the rod analyzed at
% the x-values at the next timestep. X and U_sol are set up in the same way
% and are the calculated values for RK4 and ode45 respectively. The step
% sizes for both x and t are given by delta_x and delta_t respectively.
[t_vals_heat, x_vals_heat, X_analytical_heat, X_heat, U_sol_heat, delta_x_heat, delta_t_heat] = basic_heat(t_iter_heat, x_iter_heat);

% Wave Code
t_iter_wave = 500;
x_iter_wave = 100;
c_wave = 1;

[t_vals_wave, x_vals_wave, X_analytical_wave, X_wave, delta_x_wave, delta_t_wave] = general_wave(t_iter_wave, x_iter_wave, c_wave);

t_iter_2D = 1000;
x_iter_2D = 10;
y_iter_2D = 10;
alpha_2D = 1;
[t_vals_2D, x_vals_2D, y_vals_2D, X_analytical_2D, X_2D, delta_xy_2D, delta_t_2D] = general_heat_2D(t_iter_2D, x_iter_2D, y_iter_2D, alpha_2D);