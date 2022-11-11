function value = obj_lik(y, SNRs, PHN_est, dx_est, dy_est, h_est, P)


    % update sensing matrix Omega
    Omega = [];             %sensing matrix

    % update steerring matrix
    Ar = cell(P.N,1);
    AOA = cell(P.N,1);
    for i = 1:P.N
        BS_Pos_this = P.BS_Pos(:,i);
        [Ar_n,AOA_n] = steeringMatrix(BS_Pos_this, dx_est, dy_est, P);
        Ar{i} = Ar_n;
        AOA{i} = AOA_n;
    end
    %construct sensing matrix
    for i = 1:P.N
        SNR = SNRs(i);
        tx_power = (10^(SNR/10));
        this_blk = P.x*sqrt(tx_power)*diag(exp(1i*PHN_est(:,i)))*Ar{i};
        Omega = blkdiag(Omega,this_blk);
    end
    value = norm(y-Omega*h_est)^2;
end

