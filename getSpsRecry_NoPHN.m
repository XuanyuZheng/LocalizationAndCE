function [Ar_P,Ar,AOAs] = ...
    getSpsRecry_NoPHN(dx,dy,SNRs,P)
%Get Components of the sparse recovery model;
    Ar = cell(P.N,1);
    AOAs = cell(P.N,1);


    for i = 1:P.N
        BS_Pos_this = P.BS_Pos(:,i);
        [Ar_n,AOA_n] = steeringMatrix(BS_Pos_this, dx, dy, P);
        Ar{i} = Ar_n;
        AOAs{i} = AOA_n;
    end

    % True sensing matrix
    Ar_P = [];             %true sensing matrix
    %construct actual sensing matrix
    for i = 1:P.N
        SNR = SNRs(i);
        tx_power = (10^(SNR/10));
        this_blk = P.x*sqrt(tx_power)*Ar{i};
        Ar_P = blkdiag(Ar_P,this_blk);
    end
end