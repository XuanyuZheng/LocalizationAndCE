function [cord_fixed_est,h_est,h_est_2] = loc_SRED_SPI(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P)
% Not estimating PHN, with new likelihood, with NLoS paths
%   Detailed explanation goes here

    % initialize of parameters
    % Gamma dist parameters
    beta = 0.0001;
    alpha = 0.0001;
    noise_power = 1;
    
	% initialize parameters to be estimated
    dx_est = zeros(P.Q,1);  % offset x
    dy_est = zeros(P.Q,1);  % offset y
    dphi = zeros(P.M,P.N);
    % LoS sparse channles
    plus_minus = unidrnd(2,[2*P.N*P.Q,1]);
    for i = 1:2*P.N*P.Q
       if  plus_minus(i) == 2
           plus_minus(i) = -1;
       end
    end
%     tx_power_LoS = (10.^(SNRs_LoS/10));
%     sum_powers_LoS = sum(vec(tx_power_LoS));
%     w_tilde_est = plus_minus*sqrt(sum_powers_LoS)/sqrt(P.Q*P.N)/sqrt(2); % first half real, second half img
    w_tilde_est = plus_minus/sqrt(P.Q*P.N)/sqrt(2); % first half real, second half img
%     w_tilde_est = zeros(2*P.N*P.Q,1);
    % NLoS sparse channels
%     tx_power_NLoS = (10.^(SNRs_NLoS/10));
%     sum_powers_NLoS = sum(vec(tx_power_NLoS));
%     v_tilde_est = ones(2*P.N*P.M,1)*sqrt(sum_powers_NLoS)/sqrt(P.M*P.N)/sqrt(2);               % first half real, second half img
    v_tilde_est = ones(2*P.N*P.M,1)/sqrt(P.M*P.N)/sqrt(2);               % first half real, second half img
%     v_tilde_est = zeros(2*P.N*P.M,1);
    % precomputing
    % compute covariance matrix, etc, in rad
    y_tilde = [real(y); imag(y)];
    D_y = diag(y);
    D_y_tilde = [imag(D_y);-real(D_y)];
    % covariance matrix
    W_y = eye(2*P.N*P.Nr) + 1/sqrt(noise_power) * D_y_tilde * 2*deg2rad(PHN_deg) * D_y_tilde.';
    % W_y^(-1/2)
    W_chol_inv = chol(inv(W_y));
    y_w = W_chol_inv * y_tilde;
    
    % update Omega_tilde
    [Omega,Ars,AOAs_qs] = ...
    getSpsRecry_LoS(dx_est,dy_est,SNRs_LoS,P);
    Omega_tilde = [real(Omega) -imag(Omega); imag(Omega) real(Omega)];
    W_Omega_tilde = W_chol_inv * Omega_tilde;
    
    % update Phi_tilde
    [Phi,Brs,AOAs_phis] = ...
        getSpsRecry_NLoS(dphi,SNRs_NLoS,P);
    Phi_tilde = [real(Phi) -imag(Phi); imag(Phi) real(Phi)];
    W_Phi_tilde = W_chol_inv * Phi_tilde;

    %Algorithm specifications
    bound = 1/2 * chi2inv(0.99, 2*P.N*P.Nr);
    thres_fix = 1e-6;               % stopping criteria for h for fixed grid
    iter = 1;
    max_iter = 80;                  % max number of iteration
    err_fixed = Inf;
    
    while iter <= max_iter && err_fixed >= thres_fix
        
        w_tilde_est_old = w_tilde_est;
        v_tilde_est_old = v_tilde_est;
        u_tilde_est_old = [w_tilde_est_old;v_tilde_est_old];
        
        lambda = zeros(P.Q,1);
        for q = 1:P.Q
            w_bar_sum_over_n = 0;
            for n = 1:2*P.N
                w_bar_sum_over_n = w_bar_sum_over_n + norm(w_tilde_est(   (n-1)*P.Q+q   ),2)^2;       %Here use prior norm????????????
            end
            lambda(q) = sqrt(alpha+2*P.N)/(beta+sqrt(w_bar_sum_over_n));
        end
        
        %solve w, do variable substution, s.t. x(n,q) = lambda(q)*h(n,q)
        A = zeros(2*P.N*P.Nr,2*P.N*P.Q);
        for n1 = 1:2*P.N
           for q1 = 1:P.Q
              A(:,(n1-1)*P.Q+q1) = W_Omega_tilde(:,(n1-1)*P.Q+q1)/(lambda(q1));
           end
        end
        
        gamma = zeros(P.N*P.M,1);
        for nm = 1:P.N*P.M
            v_sum_overNM = v_tilde_est(nm).^2 + v_tilde_est(nm+P.N*P.M).^2;
            gamma(nm) = sqrt(alpha+2)/(beta+sqrt(v_sum_overNM));
        end
        
        % solve v, do variable substution, s.t. z(n,m) = gamma(n,m)*v(n,m)
        B = zeros(2*P.N*P.Nr,2*P.N*P.M);
        for n2 = 1:2
           for mn = 1:P.M*P.N
              B(:,(n2-1)*P.M*P.N+mn) = W_Phi_tilde(:,(n2-1)*P.M*P.N+mn)/(gamma(mn));
           end
        end
        
        cvx_begin quiet
%                 cvx_solver mosek


            variable x1(2*P.N,P.Q) 
            variable x2(2,P.N*P.M) 
            
            minimize( sum(norms(x1)) +  sum(norms(x2))  )
            subject to
            sum_square_abs(y_w-A*vec(transpose(x1)) - B*vec(transpose(x2))) <= (1*bound); %#ok<VUNUS>

%             variable x1(2*P.N,P.Q) complex
%             minimize( sum(norms(x1)) + sum_square_abs(y1-A*vec(transpose(x1)))   )
        cvx_end
        
        % recover w
        for n1 = 1:2*P.N
           for q1 = 1:P.Q
              w_tilde_est((n1-1)*P.Q+q1) = x1(n1,q1)/(lambda(q1));
           end
        end
        
        % recover v
        for n2 = 1:2
           for mn = 1:P.M*P.N
                v_tilde_est((n2-1)*P.M*P.N+mn) = x2(n2,mn)/gamma(mn);
           end
        end
        
        u_tilde_est = [w_tilde_est;v_tilde_est];
        err_fixed = norm(u_tilde_est-u_tilde_est_old)/norm(u_tilde_est_old)
        iter = iter + 1;
        
    end
    
    w_est = w_tilde_est(1:P.N*P.Q) + 1i * w_tilde_est(P.N*P.Q+1:end);
    w_est_mat = zeros(P.Q,P.N);
    for i = 1:P.N
       for j = 1:P.Q
           w_est_mat(j,i) = w_est((i-1)*P.Q+j);
       end
    end
    w_est_mat_sum = sum(abs(w_est_mat'));
    figure()
    stem((w_est_mat_sum'))
    position_est = find(w_est_mat_sum==max(w_est_mat_sum));
    cord_fixed_est = P.qs{position_est};
    
    v_est = v_tilde_est(1:P.N*P.M) + 1i * v_tilde_est(P.N*P.M+1:end);
    figure()
    stem(abs(v_est))
    
    h_NLoS_est = reshape(Phi * v_est,[P.Nr,P.N]);
    h_LoS_est = reshape(Omega * w_est,[P.Nr,P.N]);
    h_est = h_LoS_est + h_NLoS_est;
    
    supp_w = position_est : P.Q : position_est+(P.N-1)*P.Q;
    w_est_1_supp = zeros(size(w_est));
    w_est_1_supp(supp_w) = w_est(supp_w);
    h_LoS_est_2 = reshape(Omega * w_est_1_supp,[P.Nr,P.N]);
    h_est_2 = h_LoS_est_2 + h_NLoS_est;
    


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end

