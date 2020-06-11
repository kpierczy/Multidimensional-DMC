%==============================================================
% Returns object's output for a given input vector. 
%
% @param U : input vector
% @returns : object's output
%
% @note : Function uses row, not column vectors.
%==============================================================
function Y = simulate(U)

    % Get length of the simulation 
    simulation_time = size(U,2);

    % Initialize output vector
    Y = zeros(1, simulation_time);

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

