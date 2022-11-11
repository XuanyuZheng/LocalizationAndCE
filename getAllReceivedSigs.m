function [y,PHN,w,v,h_tr,AoA_NLoS_tr,CovMat] = ...
    getAllReceivedSigs(MU_Pos,scat_pos,SNRs_LoS,SNRs_NLoS,PHN_deg,P)
%getAllReceivedSigs
% Input MU position, scatterers positions, receive SNRs, PHN_std and P
% struct
% Output: received signals across all BSs, true PHNs, N LoS channel
% coefficients, N x K NLoS channel coefficients, reconstructed LOS+NLoS overall channel
    K = size(scat_pos, 2);                          % scat_num
    y = zeros(P.N*P.Nr,1);                          % stacked received signal of the BS in the sequence of BS_Pos, C^(Nr,N)
    PHN = zeros(P.Nr,P.N);                          % real PHN, record
    w = zeros(P.N,1);                               % N LoS channel gains
    v = zeros(P.N,K);                               % N x K NLoS channel gains
    h_tr = zeros(P.Nr,P.N);
    AoA_NLoS_tr = zeros(P.N,K);
%     q = MU_Pos;
    for n = 1:P.N
        SNR_LoS = SNRs_LoS(n);                      % LoS SNR to n-th BS
        SNR_NLoS = SNRs_NLoS(n,:);                  % NLoS SNRs to n-th BS
        p = P.BS_Pos(:,n);
%         [yn, PHN_n, h_n, CovMat] = getReceivedSigBS(p, MU_Pos, scat_pos, SNR_LoS, SNR_NLoS, PHN_deg, P);
        [yn, PHN_n, w_n, v_n, h_tr_n, ~, AoA_NLoS_n, CovMat] = getReceivedSigBS(p,MU_Pos,scat_pos,SNR_LoS,SNR_NLoS,PHN_deg,P);
        w(n) = w_n;
        v(n,:) = v_n;
        PHN(:,n) = PHN_n;
        AoA_NLoS_tr(n,:) = AoA_NLoS_n;
        y((n-1)*P.Nr + 1 : n*P.Nr) = yn;
        h_tr(:,n) = h_tr_n;
    end
end

