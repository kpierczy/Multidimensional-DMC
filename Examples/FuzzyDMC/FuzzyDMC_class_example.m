%========================================================================
% File name   : FuzzyDMC_class_example.m
% Date        : 11th June 2020
% Author      : Krzysztof Pierczyk
% Description : Simple example that shows how to use FuzzyDMC class
%========================================================================

clear
clc

% First argument passed to the FuzzyDMC constructor is a cell array
% of structures containing regulator-specific parameters like horizons,
% or lambda values. In this example we will create 4 regulater with the
% same parameters that will differ in membership functions that
% describes them. Shape of the example object is [1, 0, 1] (SISO)
%
% @note : rss is an abbreviation from 'regulator specific structure'

% Set number of regulators
regulators_num = 10;

% Determine DMCs parameters (in this example the same for all regulators)
params_t.N      = 14;
params_t.Nu     = 14;
params_t.D      = 50;
params_t.Dz     = 0;
lambda          = 1;
psi             = 1;

% @note : lambda and psi matrices are not contained in the parameters
%         struct. It is because they have to be set after initializaing
%         internal object's model of the local regulator. To get more info 
%         @see DMC doc.

% Stack parameters structure into the cell array (in general parameters
% of local regulators could be different)
params = cell(regulators_num, 1);
for i = 1:regulators_num
   params{i} = params_t; 
end

% Next step is to choose shape of the fuzzy sets used. This parameter
% is common for all local regulators. Parameters of the same sets
% are individual for each local regulator and will be set in the
% underlying loop
membership_fun = 'gaussian';

% note : in this case we use default 'gaussian' sets sizes
%        @see FuzzyDMC to check whether other shapes are
%        available

% Create cell array of 10 parameters sets for local models of the object
loaclModels = cell(regulators_num, 1);
% For each regulator we define different parameters concerning
% process of step responses gathering and fuzzy sets shapes
for i = 1:regulators_num
    
    % Models can be defined given a function which meets some requirements
    % (@see DMC.updateModelFunction()) or with step responses gathered from
    % the object (@see DMC.updateModel()). In this case we use function
    loaclModels{i}.object = @(U)(simulate(U));
    
    % Shape of the model is vector with number of CVs, DVs and PVs (in
    % this order)
    loaclModels{i}.shape = [1; 0; 1];
    
    % Each subregulator will gather step response from the object to
    % prepare own local object's model. In this step we set work point
    % that CVs steps will be performed (in this case we have only one
    % CV so vector converges to scalar)
    loaclModels{i}.init_point_u = -1 + (i-1)*2/(regulators_num - 1);
    
    % Size of the step sizes used during step response gathering is
    % also configurable
    if i == regulators_num
        loaclModels{i}.step_size_u = -0.2;
    else
        loaclModels{i}.step_size_u = 0.2;
    end    
   
    % At last, we need to configure parameters of the fuzzy sets used to
    % fuzzify local models' CVs. The msot reasonable tactic seems to be
    % setting center of the fuzzy set at the work point that step response
    % will be gathered from (or maybe an average of the work point
    % and work point + step size (?))
    loaclModels{i}.c     = loaclModels{i}.init_point_u;
    
    % Sigma is parameter that decides how spreaded fuzzy sets will be
    loaclModels{i}.sigma = 0.15;    
    
end

% Now we can construct out regulator
FuzzyDMC_reg = FuzzyDMC(params, loaclModels, membership_fun);

% Initializa 'lambda' and 'psi' matrices
for i = 1:regulators_num
   FuzzyDMC_reg.regulators{i}.lambda = lambda;
   FuzzyDMC_reg.regulators{i}.psi    = psi;
end

% @note : In this example all local regulator's 'lambda' and 'psi'
%         matrices are set to (default) scalar value 1. If 2D matrices
%         are used, their sizes should be adjusted to the particular
%         regulator's parameters AFTER regulators' creation. To get
%         know why it is so @see DMC doc.

%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

% Simulation length
sim_time = 700;

% Lets initialize some desired trajectory
y_zad_trajectory_struct.size  = sim_time;
y_zad_trajectory_struct.steps = [1, -0.25; 100, -0.5; 200, -1; 300, -2; ...
                                 400, -0.25; 500, -2.5; 600, -1];
y_zad = make_trajectory(y_zad_trajectory_struct);

% Initialize PVs vector
y = zeros(sim_time, 1);

% Initialize CVs vector
u = zeros(sim_time, 1);

for i=1:sim_time
    y(i) = simulate_step(u, y, i);
    u(i) = FuzzyDMC_reg.compute(y_zad(i:end), y(i));
end


%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

hold on
plot(y)
plot(y_zad)
hold off

figure
hold on
ls = linspace(-1, 1, 200);
for i = 1:regulators_num
    plot(gaussmf(ls, [FuzzyDMC_reg.sigma(i), FuzzyDMC_reg.c(i)]))
end

%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

clearvars -except FuzzyDMC_reg y  y_zad u
