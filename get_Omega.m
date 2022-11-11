function Omega = get_Omega(Ar,PHN,SNRs,P)
    % True sensing matrix
    Omega = [];             %true sensing matrix
    %construct actual sensing matrix
    for i = 1:P.N
        SNR = SNRs(i);
        tx_power = (10^(SNR/10));
        this_blk = P.x*sqrt(tx_power)*diag(exp(1i*PHN(:,i)))*Ar{i};
        Omega = blkdiag(Omega,this_blk);
    end
end

