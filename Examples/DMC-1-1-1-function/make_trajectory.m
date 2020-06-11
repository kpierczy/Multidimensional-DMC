%========================================================================
% File name   : DMC_class_example.m
% Date        : 27th April 2020
% Author      : Krzysztof Pierczyk
% Description : Make trajectory function returns vector of trajectory 
%               values constructed with values given in a 
%               'trajectory_struct' structure.
%
% @param trajectory_struct : structure describing the trajectory
%        with following fields:
%
%        .size - number of elements in the trajectory vector
%        
%        .steps - matrix (M x 2) containing trajectory values
%                 steps. The first column contains elements
%                 when the steps occurs, when the second holds
%                 values of the trajectory after steps.
%
% @note : Trajectory vector is initialized with zeros, so if the
%         first step is not applied at the 0th element, first
%         (N - 1) elements are set to 0 (where N is values hold in
%         trjaectory_struct.steps(1, 1))
%========================================================================

function [trajectory] = make_trajectory(trajectory_struct)

    % Check if structure contains appropriate fields
    if not(isfield(trajectory_struct, 'size'))
       error('Structure does not contains trajectory size!') 
    end
    if not(isfield(trajectory_struct, 'steps'))
       error('Structure does not contains trajectory values!') 
    end

    % Check size of the trajectory
    if trajectory_struct.size < 1
       error('Trajectory size should be greater than 0!') 
    end
    
    % Check number of trajectory value steps
    if size(trajectory_struct.steps, 1) > trajectory_struct.size
       error('Too many steps!') 
    end

    % Initialize trajectory vector
    trajectory = zeros(trajectory_struct.size, 1);
    
    % Set appropriate step values
    if size(trajectory_struct.steps, 1) ~= 0
        for i=1:size(trajectory_struct.steps, 1)
            trajectory(trajectory_struct.steps(i, 1):size(trajectory, 1)) = trajectory_struct.steps(i, 2);
        end
    end

end

