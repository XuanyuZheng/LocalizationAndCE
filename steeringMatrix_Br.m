function [Br,AOAs] = steeringMatrix_Br(dphi_n, P)
%   steeringMatrix_Br  Generate NLoS steering matrix Br for a particular BS
%   dphi_n is the off-grid for phi at BS n, a M x 1 vector

    % construct the coordinates of the (grids + offsets)
    q_phi = zeros(P.M,1);      % coordinates of the fixed grids + offset

    for m = 1:P.M
        q_phi(m) = P.phi_mesh(m) + dphi_n(m);
    end
    
    Br = zeros(P.Nr,P.M);
    antenna_pos = P.antenna_pos;
    AOAs = zeros(P.M,1);
    for m = 1:P.M
        theta = q_phi(m);
        % construct the steering vector
        a = steering_vec_NLoS(theta,antenna_pos,P.Nr,P.lambda_c);
        Br(:,m) = a;
        AOAs(m) = theta;
    end

end
