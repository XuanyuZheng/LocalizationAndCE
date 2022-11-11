function [Omega,Ar,AOAs] = ...
    getSpsRecry_LoS(dx,dy,SNRs,P)
%Get Components of the sparse recovery model;
    Ar = cell(P.N,1);
    AOAs = cell(P.N,1);


    for i = 1:P.N
        BS_Pos_this = P.BS_Pos(:,i);
        [Ar_n,AOA_n] = steeringMatrix_Ar(BS_Pos_this, dx, dy, P);
        Ar{i} = Ar_n;
        AOAs{i} = AOA_n;
    end

    % True sensing matrix
    Omega = [];             %true sensing matrix
    %construct actual sensing matrix
    for i = 1:P.N
        SNR = SNRs(i);
        tx_power = (10^(SNR/10));
        this_blk = P.x*sqrt(tx_power)*Ar{i};
%         this_blk = P.x*Ar{i};
        Omega = blkdiag(Omega,this_blk);
    end
end

