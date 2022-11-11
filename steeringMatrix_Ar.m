function [Ar,AOAs] = steeringMatrix_Ar(BS_Pos_this, dx, dy, P)
%steeringMatrix Generate steering matrix for a particular BS
%   BS_Pos_this is the position of the BS
%   the total_length x total_length square is meshed to grid_side x grid_side small squares

    % construct the coordinates of the (grids + offsets)
    q = cell(P.grid_side,P.grid_side);      % coordinates of the fixed grids + offset
    for i = 1:P.grid_side
       for j = 1:P.grid_side
           q{i,j} = P.qs_mesh{i,j} + [dx((i-1)*P.grid_side+j);dy((i-1)*P.grid_side+j)];
       end
    end
    q = vec(q');
    
    Ar = zeros(P.Nr,P.Q);
    antenna_pos = P.antenna_pos;
    AOAs = zeros(P.Q,1);
    for j = 1:P.Q
        % construct the steering vector
        [a,theta] = steering_vec(BS_Pos_this,q{j},antenna_pos,P.Nr,P.lambda_c);
        Ar(:,j) = a;
        AOAs(j) = theta;
    end

end

