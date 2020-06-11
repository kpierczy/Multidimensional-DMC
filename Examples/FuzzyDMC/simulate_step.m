%==============================================================
% Returns object's output for a given input vector for
% a required discrete moment. If sizes of the input or output
% vector is lower than the required moment the error is thrown.
%
% @param U : input vector
% @param Y : output vector
% @param moment : required moment
% @returns : object's output at the 'moment'th momenht 
%
%==============================================================
function y =simulate_step(U, Y, moment)

    % Throw error when vectors' sizes are not equal
    if moment > size(U,1) || moment > size(Y,1)
        error('Required moment exceeds given vectors sizes!')
    end

    % Simulate
    if ( moment > 6)
        y = object(U(moment - 5), U(moment - 6), Y(moment - 1), Y(moment - 2));
    elseif ( moment > 5 )
        y = object(U(moment - 5),             0, Y(moment - 1), Y(moment - 2));
    elseif ( moment > 2 )
        y = object(            0,             0, Y(moment - 1), Y(moment - 2));
    elseif ( moment > 1 )
        y = object(            0,             0, Y(moment - 1),             0);
    else
        y = object(            0,             0,             0,             0);
    end

end

