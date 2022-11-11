function [Phi,Br,AOAs] = ...
    getSpsRecry_NLoS(dphi,SNRs,P)
%Get Components of the sparse recovery model;
    Br = cell(P.N,1);
    AOAs = cell(P.N,1);


    for n = 1:P.N
        dphi_n = dphi(:,n);
        [Br_n,AOAs_phis_n] = steeringMatrix_Br(dphi_n, P);
        Br{n} = Br_n;
        AOAs{n} = AOAs_phis_n;
    end

    % True sensing matrix
    Phi = [];             %true sensing matrix
    %construct actual sensing matrix
    for n = 1:P.N
        SNR = SNRs(n,1);
        tx_power = (10^(SNR/10));
        this_blk = P.x*sqrt(tx_power)*Br{n};
%         this_blk = P.x*Br{n};
        Phi = blkdiag(Phi,this_blk);
    end
end

