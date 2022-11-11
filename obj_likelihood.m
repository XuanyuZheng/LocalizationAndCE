function value = obj_likelihood(y_tilde, SNRs_LoS,SNRs_NLoS, ...
    dx_est, dy_est, dphi_est, w_tilde_est, v_tilde_est, W_chol_inv, P)


    % update Omega_tilde
    [Omega,~,~] = ...
    getSpsRecry_LoS(dx_est,dy_est,SNRs_LoS,P);
    Omega_tilde = [real(Omega) -imag(Omega); imag(Omega) real(Omega)];
    
    % update Phi_tilde
    [Phi,~,~] = ...
        getSpsRecry_NLoS(dphi_est,SNRs_NLoS,P);
    Phi_tilde = [real(Phi) -imag(Phi); imag(Phi) real(Phi)];
    
    value = norm(W_chol_inv*y_tilde - ...
        W_chol_inv*Omega_tilde*w_tilde_est - W_chol_inv*Phi_tilde*v_tilde_est)^2;

    

%     value = norm(y-Omega*h_tilde_est)^2;
    
end