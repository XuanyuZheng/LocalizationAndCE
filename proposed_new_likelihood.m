function [cord_fixed_est, h_est_N, AOAs_est] = proposed_new_likelihood(y,SNRs,PHN_deg,PHN_CovMat,P)
% Not estimating PHN, with new likelihood

    % Gamma dist parameters
    beta = 0.0001;
    alpha = 0.01;
    noise_power = 1;
    
	% initialize parameters to be estimated
    dx_est = zeros(P.Q,1);  % offset x
    dy_est = zeros(P.Q,1);  % offset y
%     PHN = zeros(P.Nr,P.N);   % fixed to zero
    
    plus_minus = unidrnd(2,2*P.N*P.Q,1);
    for i = 1:2*P.N*P.Q
       if  plus_minus(i) == 2
           plus_minus(i) = -1;
       end
    end
    h_tilde_est = plus_minus/sqrt(P.Q)/sqrt(2); % first half real, second half img
%     h_tilde_est = plus_minus*0;
    %Algorithm specifications
    thres_fix = 1e-4;               % stopping criteria for h for fixed grid
    err_fix = Inf;                      % First error
    iter = 1;
    max_iter = 30;                  % max number of iteration
    
    y_tilde = [real(y); imag(y)];
    D_y = diag(y);
    D_y_tilde = [imag(D_y);-real(D_y)];
%     W_y = eye(2*P.N*P.Nr) + 1/sqrt(noise_power) * D_y_tilde * 2*deg2rad(PHN_deg) * D_y_tilde.';
    W_y = eye(2*P.N*P.Nr) + 1/sqrt(noise_power) * D_y_tilde * 2*PHN_CovMat * D_y_tilde.';
    W_chol_inv = chol(inv(W_y));
    y1 = W_chol_inv * y_tilde;
    
    % update Omega_tilde
    
    [Ar_P,~,AOAs] = ...
    getSpsRecry_NoPHN(dx_est,dy_est,SNRs,P);
    Omega_tilde = [real(Ar_P) -imag(Ar_P); imag(Ar_P) real(Ar_P)];
    W_Omega_tilde = W_chol_inv * Omega_tilde;
    
    

    bound = 1/2 * chi2inv(0.99, 2*P.N*P.Nr);
    
    
    while iter <= max_iter && err_fix >= thres_fix
        h_tilde_est_old = h_tilde_est;
       
        % update weights
        lambda = zeros(P.Q,1);
        for q = 1:P.Q
            h_bar_sum_over_n = 0;
            for n = 1:2*P.N
                h_bar_sum_over_n = h_bar_sum_over_n + norm(h_tilde_est(   (n-1)*P.Q+q   ),2)^2;       %Here use prior norm????????????
            end
            lambda(q) = (alpha+P.N)/(beta+sqrt(h_bar_sum_over_n));
        end
        %solve h
        % do variable substution, s.t. x(n,q) = lambda(q)*h(n,q)
        A = zeros(2*P.N*P.Nr,2*P.N*P.Q);
        for n1 = 1:2*P.N
           for q1 = 1:P.Q
              A(:,(n1-1)*P.Q+q1) = W_Omega_tilde(:,(n1-1)*P.Q+q1)/(lambda(q1));
           end
        end
        
        cvx_begin quiet
            cvx_solver mosek
            
            
            variable x1(2*P.N,P.Q) complex
            minimize( sum(norms(x1))   )
            subject to
            sum_square_abs(y1-A*vec(transpose(x1))) <= (1*bound); %#ok<VUNUS>
            
%             variable x1(2*P.N,P.Q) complex
%             minimize( sum(norms(x1)) + sum_square_abs(y1-A*vec(transpose(x1)))   )
            
            
        cvx_end

        for n1 = 1:2*P.N
           for q1 = 1:P.Q
              h_tilde_est((n1-1)*P.Q+q1) = x1(n1,q1)/(lambda(q1));
           end
        end
        
        err_fix = norm(h_tilde_est-h_tilde_est_old)/norm(h_tilde_est_old);
        iter = iter + 1;
        
        h_est = h_tilde_est(1:P.N*P.Q) + 1i * h_tilde_est(P.N*P.Q+1:end);
        h_est_mat = zeros(P.Q,P.N);
        for i = 1:P.N
           for j = 1:P.Q
               h_est_mat(j,i) = h_est((i-1)*P.Q+j);
           end
        end
        h_est_mat_sum = sum(abs(h_est_mat'));
%         figure()
%         plot((h_est_mat_sum'))
        position_est = find(h_est_mat_sum==max(h_est_mat_sum));          % estimated position, recorded!
%         fprintf('iter=%d, position_est = %d, dx_est = %f, dy_est = %f, error_fixed = %f\n',...
%             iter,position_est,dx_est(position_est),dy_est(position_est),err_fix);
    end
    
    
    h_est = h_tilde_est(1:P.N*P.Q) + 1i * h_tilde_est(P.N*P.Q+1:end);
    h_est_mat = zeros(P.Q,P.N);
    for i = 1:P.N
       for j = 1:P.Q
           h_est_mat(j,i) = h_est((i-1)*P.Q+j);
       end
    end
    h_est_mat_sum = sum(abs(h_est_mat'));
%     figure()
%     stem((h_est_mat_sum'))
    position_est = find(h_est_mat_sum==max(h_est_mat_sum));          % estimated position, recorded!

    cord_fixed_est = [(mod(position_est,P.grid_side)-1)*P.grid_length+P.grid_length/2-P.grid_side*P.grid_length/2;...
    -((ceil(position_est/P.grid_side)-1)*P.grid_length+P.grid_length/2)+P.grid_side*P.grid_length/2];

    if mod(position_est,P.grid_side) == 0         %deal with mod = 0 but 
        cord_fixed_est(1) = (P.grid_side-1)*P.grid_length+P.grid_length/2-P.grid_side*P.grid_length/2;
    end
    h_est_N = h_est_mat(position_est, :).';
    
    supp = position_est : P.Q : position_est+(P.N-1)*P.Q;
    Atmp = cell2mat(AOAs);
    AOAs_est = Atmp(supp);
%     fprintf('\ncord_fixed_est = [%f, %f]\n',cord_fixed_est(1), cord_fixed_est(2));
end

