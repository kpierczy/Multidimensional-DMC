%========================================================================
% File name   : FuzzyDMC_class_example.m
% Date        : 28th April 2020
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
rss.N      = 14;
rss.Nu     = 14;
rss.D      = 50;
rss.Dz     = 0;
rss.lambda = 1;

% Create cell array of 4 identical regulators
regulators = cell(regulators_num, 1);

% For each regulator we define different parameters concerning
% process of step responses gathering and fuzzy sets shapes
for i = 1:regulators_num
    regulators{i, 1} = rss;
    regulators{i, 1}.sr_point = [-1 + (i-1)*2/(regulators_num - 1)];
    if i == regulators_num
        regulators{i, 1}.sr_step_size = [-0.2];
    else
        regulators{i, 1}.sr_step_size = [0.2];
    end    
    regulators{i, 1}.c     = -1 + (i-1)*2/(regulators_num - 1);
    regulators{i, 1}.sigma = 0.15;    
end

% The second parameter of the constructor is a structure containing
% informations that are common to all regulators
%
% @note : cs is and abbreviation from 'common structure'

% Shape of the object is vector containing number of CVs, Dvs and PVs
cs.shape = [1, 0, 1];

% Handle to the function representing object (to see function's
% requirements @see DMC(...) constructor)
cs.object = @(U)(simulate(U));

% Type of the membership function - athis time only 'gaussian' is available
cs.membership_fun = 'gaussian';

% Now we can construct out regulator
FuzzyDMC_reg = FuzzyDMC(regulators, cs);

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

% figure
% hold on
% ls = linspace(-1, 1, 200);
% for i = 1:regulators_num
%     plot(gaussmf(ls, [FuzzyDMC_reg.sigma(i), FuzzyDMC_reg.c(i)]))
% end


%========================================================================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Cleaning %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%========================================================================

clear
