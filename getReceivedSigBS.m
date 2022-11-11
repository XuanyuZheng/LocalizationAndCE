function [y, PHN, w_n, v_n, h_tr, AoA_LoS, AoA_NLoS, CovMat] = ...
    getReceivedSigBS(p,q,scat_pos,SNR_LoS,SNR_NLoS,PHN_deg,P)
    %RecAtBS: This function generates the received signal at one BS
    % P.x: unit transmitted symbol, C^(1,1)
    % p: (px,py) position of this BS, R^(2,1)
    % q: (qx,qy) position of the MU, R^(2,1)
    % scat_pos: position of the scatteres R^(2,K)
    % Nr: The number of antennas on each BS, R^(1,1)
    % lambda_c: carrier frequency, R^(1,1)
    % SNR_LoS: in terms of power, R^(1,1)
    % SNR_NLoS: in terms of power, R^(1,K)
    % Lambda: PHN matrix of this BS
    K = size(scat_pos,2);       % number of scatterers
    
    % noise and multipath powers
    noise_power = 1;
    noise_std = sqrt(noise_power);
    noise = noise_std * sqrt(1/2) * (randn(P.Nr,1) + 1i * randn(P.Nr,1));     % noise variance is noise_power
    
    tx_power_LoS = (10^(SNR_LoS/10));
    tx_power_NLoS = (10.^(SNR_NLoS/10));
    
    w_n = sqrt(1/2)*(randn + 1i*randn);                % LoS channel gain realization of this basestation
    v_n = sqrt(1/2)*(randn(1,K) + 1i*randn(1,K));      % NLoS channel gains realization of this basestation

    std_dev = deg2rad(PHN_deg);
    PHN = std_dev*randn(P.Nr,1);                      %Gaussian dis
%     PHN = -std_dev + (2*std_dev)*rand(P.Nr,1);          %Uniform dis
    Lambda = diag(exp(1i*PHN));
    
    len = P.Nr;
    cor = 0.2;
    CovMat = zeros(len)*cor;
    CovMat(logical(eye(len))) = ones(len,1);
%     L = chol(CovMat);

%     s = randn(len,1);
%     PHN = L' * s * std_dev;
%     Lambda = diag(exp(1i*PHN));
%     Cov_re = cov(r.');

    % construct the steering vector for LoS
    [a_LoS, AoA_LoS] = steering_vec(p,q,P.antenna_pos,P.Nr,P.lambda_c);
    
    
    
    % construct the steering vector for NLoS
    AoA_NLoS = zeros(1,K);
    a_NLoS = zeros(P.Nr,K);                                   %steering vector
    for k = 1:K
        q_sc = scat_pos(:,k);                                     % coordinate of the k-th scatterer
        [a_NLoS_this, AOA_NLoS_this] = steering_vec(p,q_sc,P.antenna_pos,P.Nr,P.lambda_c);
        a_NLoS(:,k) = a_NLoS_this;
        AoA_NLoS(k) = AOA_NLoS_this;
    end
    
    h_LoS_tr = sqrt(tx_power_LoS)*w_n*a_LoS;
    h_NLoS_tr = zeros(P.Nr,1);
    for k = 1:K
        h_NLoS_tr = h_NLoS_tr + sqrt(tx_power_NLoS(k))*v_n(k)*a_NLoS(:,k);
    end
    h_tr = h_LoS_tr + h_NLoS_tr;
    
    y =  Lambda*(h_tr)*P.x + noise;

end