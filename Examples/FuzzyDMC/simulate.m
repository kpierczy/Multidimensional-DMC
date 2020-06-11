%==============================================================
% Returns object's output for a given input vector. 
%
% @param U : input vector of size (X x 1)
% @returns : object's output of size (X x 1)
%==============================================================
function Y = simulate(U)

    % Get length of the simulation 
    simulation_time = size(U,1);

    % Initialize output vector
    Y = zeros(simulation_time, 1);

    % Simulate
    for i = 1:simulation_time

        if (i > 6)
            Y(i) = object(U(i - 5), U(i - 6), Y(i - 1), Y(i - 2));
        elseif ( i > 5 )
            Y(i) = object(U(i - 5),        0, Y(i - 1), Y(i - 2));
        elseif ( i > 2 )
            Y(i) = object(       0,        0, Y(i - 1), Y(i - 2));
        elseif ( i > 1 )
            Y(i) = object(       0,        0, Y(i - 1),        0);
        else
            Y(i) = object(       0,        0,        0,        0);
        end

    end
end

