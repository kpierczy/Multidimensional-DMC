%========================================================================
% File name   : DMC_class_example.m
% Date        : 27th April 2020
% Author      : Krzysztof Pierczyk
% Description : Simple example that shows how to use DMC class
%========================================================================

clear
clc

% First argument passed to the DMC constructore is structure that
% describes how step responses of the CV-PV and DV-PV tracks
% should be gathered to constructthe regulator. We need to initializa
% it with certain fields

% Handle to the function representing object (to see function's
% requirements @see DMC(...) constructor)
step_responses_struct.object = @(U, Z)(simulate(U, Z));

% Shape of the object is vector with number of [CVs, DVs, PVs]
step_responses_struct.shape = [1, 1, 1];

% Values of CVs and DVs used to gather step responses
step_responses_struct.init_point_u = [0];
step_responses_struct.init_point_z = [0];

% Sizes of CVs and DVs steps used to gather step responses
step_responses_struct.step_size_u = [0.2];
step_responses_struct.step_size_z = [0.2];

% Tolerance - @see DMC(...) constructor
step_responses_struct.tol = 0.0001;


% The second argument passed to the constructor is structure containing
% DMC's parameters. Parameters should meet some natural requirements
% (e.g. be integers, be positive, non-negative, etc.). If horizons
% and dynamic ranks are set to high, DMC constructor will cut them
% to the maximum available values.

param_struct.N      = 100;
param_struct.Nu     = 30;
param_struct.D      = 200;
param_struct.Dz     = 39;
param_struct.lambda = 0.5;

% Now we can create our regulator
DMC_reg = DMC(step_responses_struct, param_struct);

%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

% Simulation length
sim_time = 400;

% Lets initialize some desired trajectory
y_zad_trajectory_struct.size  = sim_time;
y_zad_trajectory_struct.steps = [10, 1.0; 110, 3; 210, 1.5; 310, 0];
y_zad = make_trajectory(y_zad_trajectory_struct);

% Lets initialize some desired trajectory
z_trajectory_struct.size  = sim_time;
z_trajectory_struct.steps = [1, 1; 60, 0; 160, 1; 260, 0; 360, 1];
z = make_trajectory(z_trajectory_struct);

% Initialize PVs vector
y = zeros(sim_time, 1);

% Initialize CVs vector
u = zeros(sim_time, 1);

for i=1:sim_time
    y(i) = simulate_step(u, z, y, i);
    u(i) = DMC_reg.compute(y_zad(i:end), y(i), z(i));
end


%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

figure
hold on
plot(y)
plot(y_zad)


%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

clear
