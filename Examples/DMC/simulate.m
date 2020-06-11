% Returns object's output for a given input and noise
% vector. If sizes of the input and noise vectors are
% not the same, the error is thrown.
%
% @param U : input vector
% @param Z : noise vector
% @returns : object's output
%
function Y = simulate(U, Z)

% Throw error when vectors' sizes are not equal
if size(U,2) ~= size(Z,2)
    error('Lengths of both trajectories are not equal!')
end
simulation_time = size(U,2);

% Initialize output vector
Y = zeros(1, simulation_time);

% Simulate
for i = 1:simulation_time
    
    if (i > 8)
        Y(i) = object(U(i-7), U(i-8), Z(i-2), Z(i-3), Y(i-1), Y(i-2));
    elseif ( i > 7 )
        Y(i) = object(U(i-7), 0, Z(i-2), Z(i-3), Y(i-1), Y(i-2));
    elseif ( i > 3 )
        Y(i) = object(0, 0, Z(i-2), Z(i-3), Y(i-1), Y(i-2));
    elseif ( i > 2 )
        Y(i) = object(0, 0, Z(i-2), 0, Y(i-1), Y(i-2));
    elseif ( i > 1 )
        Y(i) = object(0, 0, 0, 0, Y(i-1), 0);
    else
        Y(i) = object(0, 0, 0, 0, 0, 0);
    end

end

end

