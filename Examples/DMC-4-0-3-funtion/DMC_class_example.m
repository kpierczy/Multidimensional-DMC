%========================================================================
% File name   : DMC_class_example.m
% Date        : 11th June 2020
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
objectStruct.object = @(U)(simulate(U));

% Shape of the object is vector with number of [CVs, DVs, PVs]
objectStruct.shape = [4, 0, 3];

% Values of CVs and DVs used to gather step responses
objectStruct.init_point_u = [0; 0; 0; 0];

% Sizes of CVs and DVs steps used to gather step responses
objectStruct.step_size_u = [0.2; 0.2; 0.2; 0.2];

% Tolerance - @see DMC(...) constructor
objectStruct.tol = 0.0001;


% The second argument passed to the constructor is structure containing
% DMC's parameters. Parameters should meet some natural requirements
% (e.g. be integers, be positive, non-negative, etc.). If horizons
% and dynamic ranks are set to high, DMC constructor will cut them
% to the maximum available values.

param_struct.N      = 107;
param_struct.Nu     = 54;
param_struct.D      = 204;
param_struct.Dz     = 0;
lambda              = eye(param_struct.Nu * 4) * 0.5;
psi                 = eye(param_struct.N  * 3) * 1;

% @note : lambda and psi matrices are not contained in the parameters
%         struct. It is because they have to be set after initializaing
%         internal object's model of the regulator. To get more info 
%         @see DMC doc.

% Now we can create our regulator
DMC_reg = DMC(param_struct);
DMC_reg.updateModelFunction(objectStruct);
DMC_reg.lambda = lambda;
DMC_reg.psi    = psi;

%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

% Simulation length
sim_time = 1600;

% Lets initialize some desired trajectory
% Lets initialize some desired trajectory
y1_zad_trajectory_struct.size  = sim_time;
y1_zad_trajectory_struct.steps = [1, 1.0; 400, 3.8; 800, 0.2; 1200, 1];
y1_zad = make_trajectory(y1_zad_trajectory_struct);

y2_zad_trajectory_struct.size  = sim_time;
y2_zad_trajectory_struct.steps = [1, 0];
y2_zad = make_trajectory(y2_zad_trajectory_struct);

y3_zad_trajectory_struct.size  = sim_time;
y3_zad_trajectory_struct.steps = [1, 0];
y3_zad = make_trajectory(y3_zad_trajectory_struct);

y_zad = [y1_zad, y2_zad, y3_zad];

% Initialize PVs vector
y = zeros(sim_time, 3);

% Initialize CVs vector
u = zeros(sim_time, 4);

for i=1:sim_time
    y(i, :) = simulate_step(u, y, i)';
    u(i, :) = DMC_reg.compute(y_zad(i, :), y(i, :))';
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

clearvars -except y y_zad u DMC_reg
