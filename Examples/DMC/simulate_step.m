% Returns object's output for a given input and noise
% vector for a required sicrete moment. If sizes of the 
% input or noise vector is lower than the required moment
% the error is thrown
%
% @param U : input vector
% @param Z : noise vector
% @param Y : output vector
% @param moment : required moment
% @returns : object's output at the 'moment'th momenht 
%
function y = simulate_step(U, Z, Y, moment)

% Throw error when vectors' sizes are not equal
if moment > size(U,1) || moment > size(Z,1) || moment > size(Y,1)
    error('Required moment exceeds given vectors sizes!')
end

if moment > 8
    y = object(U(moment - 7), U(moment - 8), ...
               Z(moment - 2), Z(moment - 3), ...
               Y(moment - 1), Y(moment - 2));
elseif moment > 7
    y = object(U(moment - 7), 0, ...
               Z(moment - 2), Z(moment - 3), ...
               Y(moment - 1), Y(moment - 2));
elseif moment  > 3 
    y = object(0, 0, ...
               Z(moment - 2), Z(moment - 3), ...
               Y(moment - 1), Y(moment - 2));
elseif moment > 2
    y = object(0, 0, ...
               Z(moment - 2), 0, ...
               Y(moment - 1), Y(moment - 2));
elseif moment > 1
    y = object(0, 0, ...
               0, 0, ...
               Y(moment - 1), 0);
else
    y = object(0, 0, ...
               0, 0, ...
               0, 0);
end

end

