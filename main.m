% These control the number of iterations in both space and time. Note that
% this is not the number of steps, which would be either of these values
% minus 1. The output matrices will be exactly these dimensions.
y_iter_heat = 5000;
x_iter_heat = 70;

% The t-values used for analysis are contained in t_vals, and the x-values
% used for analysis are contained in x_vals. X_analytical contains the
% computed values for the analytical solution. The first row contains the
% initial condition, and each subsequent row contains the rod analyzed at
% the x-values at the next timestep. X and U_sol are set up in the same way
% and are the calculated values for RK4 and ode45 respectively. The step
% sizes for both x and t are given by delta_x and delta_t respectively.
disp("Basic Heat Equation:")
[t_vals_heat, x_vals_heat, X_analytical_heat, X_heat, U_sol_heat, t_sol_heat, delta_x_heat, delta_t_heat] = basic_heat(y_iter_heat, x_iter_heat);
%%

% Wave Code
t_iter_wave = 5000;
x_iter_wave = 100;
c_wave = 0.4;
disp("General Wave Equation:")
[t_vals_wave, x_vals_wave, X_wave, delta_x_wave, delta_t_wave] = general_wave(t_iter_wave, x_iter_wave, c_wave);

% 2D Heat Code
t_iter_2D = 1000;
x_iter_2D = 10;
y_iter_2D = 10;
alpha_2D = 0.5;
disp("General 2D Heat Equation:")
[t_vals_2D, x_vals_2D, y_vals_2D, X_2D, delta_xy_2D, delta_t_2D] = general_heat_2D(t_iter_2D, x_iter_2D, y_iter_2D, alpha_2D);

% Poisson's Equation Code
x_iter_pois = 101;
y_iter_pois = x_iter_pois*2 - 1;
TOL_pois = 10^-8;
max_iter_pois = 100000;
disp("Poisson's Equation:")
[x_vals_pois, y_vals_pois, X_pois, delta_x_pois, delta_y_pois] = general_poisson(x_iter_pois, y_iter_pois, TOL_pois, max_iter_pois);

%%
% Radiator Heat Poisson's Equation Code
x_iter_pois = 51;
y_iter_pois = 101;
TOL_pois = 10^-8;
max_iter_pois = 40000;
disp("Radiator Fin Poisson's Equation:")
[x_vals_pois2, y_vals_pois2, X_pois2, delta_x_pois2, delta_y_pois2] = Radiator_fin(x_iter_pois, y_iter_pois, TOL_pois, max_iter_pois);

%%
% Radiator Heat Method of Lines

y_iter_heat = 6000;
x_iter_heat = 40;
disp("Radiator Fin MOL Equation:")
[t_vals_heat, x_vals_heat, X_heat, delta_x_heat, delta_t_heat] = Radiator_fin_MOL(y_iter_heat, x_iter_heat);
