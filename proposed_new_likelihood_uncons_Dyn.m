function [cord_est, h_est_N, AOAs_est] = proposed_new_likelihood_uncons_Dyn(y,SNRs,PHN_deg,PHN_CovMat,P)
% Not estimating PHN, with new likelihood

    % Gamma dist parameters
    beta = 0.0001;
    alpha = 0.0001;
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
    h_tilde_est = plus_minus/sqrt(P.Q*P.N)/sqrt(2); % first half real, second half img
%     h_tilde_est = plus_minus*0;
    %Algorithm specifications
    thres_fix = 1e-6;               % stopping criteria for h for fixed grid
%     err_fix = Inf;                      % First error
    iter = 1;
    max_iter = 80;                  % max number of iteration
    
    y_tilde = [real(y); imag(y)];
    D_y = diag(y);
    D_y_tilde = [imag(D_y);-real(D_y)];
    W_y = eye(2*P.N*P.Nr) + 1/sqrt(noise_power) * D_y_tilde * 2*deg2rad(PHN_deg) * D_y_tilde.';
    W_chol_inv = chol(inv(W_y));
    y1 = W_chol_inv * y_tilde;
    
    W_y_cov = eye(2*P.N*P.Nr) + 1/sqrt(noise_power) * D_y_tilde * 2*PHN_CovMat * D_y_tilde.';
    W_chol_cov_inv = chol(inv(W_y_cov));
    y_W = W_chol_cov_inv * y_tilde;
    
    % update Omega_tilde
    
    [Ar_P,Ar,AOAs] = ...
    getSpsRecry_NoPHN(dx_est,dy_est,SNRs,P);
    Omega_tilde = [real(Ar_P) -imag(Ar_P); imag(Ar_P) real(Ar_P)];
    W_Omega_tilde = W_chol_inv * Omega_tilde;
    
    

%     bound = 1/2 * chi2inv(0.99, 2*P.N*P.Nr);
    
    obj_value = obj_likelihood(y_tilde, SNRs, dx_est, dy_est, h_tilde_est, W_chol_cov_inv, P);
    err = Inf;
    fix_converged = 0;
    
    while iter <= max_iter && err >= thres_fix
        obj_value_old = obj_value;
        h_tilde_est_old = h_tilde_est;
       
        if fix_converged == 0
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
%                 cvx_solver mosek


%                 variable x1(2*P.N,P.Q) complex
%                 minimize( sum(norms(x1))   )
%                 subject to
%                 sum_square_abs(y1-A*vec(transpose(x1))) <= (1*bound); %#ok<VUNUS>

                variable x1(2*P.N,P.Q) complex
                minimize( sum(norms(x1)) + sum_square_abs(y1-A*vec(transpose(x1)))   )


            cvx_end

            for n1 = 1:2*P.N
               for q1 = 1:P.Q
                  h_tilde_est((n1-1)*P.Q+q1) = x1(n1,q1)/(lambda(q1));
               end
            end

%             [Ar_P,Ar,AOAs] = ...
%             getSpsRecry_NoPHN(dx_est,dy_est,SNRs,P);
%             Omega_tilde = [real(Ar_P) -imag(Ar_P); imag(Ar_P) real(Ar_P)];
%             W_Omega_tilde = W_chol_inv * Omega_tilde;

            err_fixed = norm(h_tilde_est-h_tilde_est_old)/norm(h_tilde_est_old);
            if err_fixed < 1e-6
                fix_converged = 1;
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
            end
            
            
        else            
            % Update offgrids
            % compute gradient for dx
            dx_grad = zeros(P.Q,1);
            for i = position_est%1:P.Q
                Ar_i_tilde = [];
                for n = 1:P.N
                    A_n_i_tilde = zeros(P.Nr,P.Q);
                    SNR = SNRs(n);
                    tx_power_n = (10^(SNR/10));

                    dtheta_ni_dx_i = (P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))/(((P.qs{i}(1)+dx_est(i))-P.BS_Pos(1,n))^2*...
                        (1+((P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))/(P.BS_Pos(1,n)-(P.qs{i}(1)+dx_est(i))))^2));
                    a_n_i = Ar{n}(:,i);
                    a_n_i_tilde = 2*pi*1i/P.lambda_c*dtheta_ni_dx_i*a_n_i.*(P.antenna_pos'*[-sin(AOAs{n}(i));cos(AOAs{n}(i))]);
                    A_n_i_tilde(:,i) = sqrt(tx_power_n)*a_n_i_tilde;
                    Ar_i_tilde = blkdiag(Ar_i_tilde,A_n_i_tilde);
                end
                Omega_i_tilde = P.x * [real(Ar_i_tilde) -imag(Ar_i_tilde); imag(Ar_i_tilde) real(Ar_i_tilde)];
                Omega_tilde = [real(Ar_P) -imag(Ar_P); imag(Ar_P) real(Ar_P)];
                dx_grad_this = -2*(y_W.')*W_chol_cov_inv*Omega_i_tilde*h_tilde_est + ...
                    2*(W_chol_cov_inv*Omega_tilde*h_tilde_est).'*(W_chol_cov_inv*Omega_i_tilde*h_tilde_est);
                dx_grad(i) = dx_grad_this;
            end

            % compute gradient for dy
            dy_grad = zeros(P.Q,1);
            for i = position_est%1:P.Q
                Ar_i_tilde = [];
                for n = 1:P.N
                    A_n_i_tilde = zeros(P.Nr,P.Q);
                    SNR = SNRs(n);
                    tx_power_n = (10^(SNR/10));
                    dtheta_ni_dy_i = 1/(((P.qs{i}(1)+dx_est(i))-P.BS_Pos(1,n))*(1+((P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))...
                        /(P.BS_Pos(1,n)-(P.qs{i}(1)+dx_est(i))))^2));
                    a_n_i = Ar{n}(:,i);
                    a_n_i_tilde = 2*pi*1i/P.lambda_c*dtheta_ni_dy_i*a_n_i.*(P.antenna_pos'*[-sin(AOAs{n}(i));cos(AOAs{n}(i))]);
                    A_n_i_tilde(:,i) = sqrt(tx_power_n)*a_n_i_tilde;
                    Ar_i_tilde = blkdiag(Ar_i_tilde,A_n_i_tilde);
                end
                Omega_i_tilde = P.x * [real(Ar_i_tilde) -imag(Ar_i_tilde); imag(Ar_i_tilde) real(Ar_i_tilde)];
                Omega_tilde = [real(Ar_P) -imag(Ar_P); imag(Ar_P) real(Ar_P)];
                dy_grad_this = -2*(y_W.')*W_chol_cov_inv*Omega_i_tilde*h_tilde_est + ...
                    2*(W_chol_cov_inv*Omega_tilde*h_tilde_est).'*(W_chol_cov_inv*Omega_i_tilde*h_tilde_est);
                dy_grad(i) = dy_grad_this;
            end

%             dx_est = dx_est - dx_grad * 10;
%             dy_est = dy_est - dy_grad * 10;

    %         value1 = obj_likelihood(y_tilde, SNRs, dx_est, dy_est, h_tilde_est, W_chol_cov_inv, P)
    %         value2 = obj_likelihood(y_tilde, SNRs, 0*dx_est, 0*dy_est, h_tilde_est, W_chol_cov_inv, P)

            % back tracking
            step_size = 100;
            alpha_bt = 0.1;
            beta_bt = 0.3;

            obj_this_y = obj_likelihood(y_tilde, SNRs, dx_est, dy_est, h_tilde_est, W_chol_cov_inv, P);
            while 1
                obj_dyn_y = obj_likelihood(y_tilde, SNRs,...
                    dx_est + step_size*(-dx_grad),...
                    dy_est + step_size*(-dy_grad),...
                    h_tilde_est, W_chol_cov_inv, P);
    %             obj_dyn_y = obj_lik(y,SNRs,PHN_est,...
    %                 dx_est + step_size*(-dx_grad),...
    %                 dy_est + step_size*(-dy_grad),h_est, P);
                if obj_dyn_y < obj_this_y + alpha_bt*step_size*[dx_grad; dy_grad]'*(-[dx_grad; dy_grad])
                    break
                else
                    step_size = beta_bt*step_size;
                end
                if step_size < 1e-30    % avoid infinite loop for wrong gradient
                   break; 
                end
            end
            dx_est = dx_est + step_size * (-dx_grad);
            dy_est = dy_est + step_size * (-dy_grad);

            [Ar_P,Ar,AOAs] = ...
            getSpsRecry_NoPHN(dx_est,dy_est,SNRs,P);
            Omega_tilde = [real(Ar_P) -imag(Ar_P); imag(Ar_P) real(Ar_P)];
            W_Omega_tilde = W_chol_inv * Omega_tilde;
            
            obj_value = obj_likelihood(y_tilde, SNRs, dx_est, dy_est, h_tilde_est, W_chol_cov_inv, P);
            err = abs(obj_value_old - obj_value)/abs(obj_value_old);
            
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
%                 cvx_solver mosek


                variable x1(2*P.N,P.Q) complex
%                 minimize( sum(norms(x1))   )
%                 subject to
%                 sum_square_abs(y1-A*vec(transpose(x1))) <= (1*bound); %#ok<VUNUS>

    %             variable x1(2*P.N,P.Q) complex
                minimize( sum(norms(x1)) + sum_square_abs(y1-A*vec(transpose(x1)))   )


            cvx_end

            for n1 = 1:2*P.N
               for q1 = 1:P.Q
                  h_tilde_est((n1-1)*P.Q+q1) = x1(n1,q1)/(lambda(q1));
               end
            end

%             [Ar_P,Ar,AOAs] = ...
%             getSpsRecry_NoPHN(dx_est,dy_est,SNRs,P);


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
%             fprintf('iter=%d, position_est = %d, dx_est = %f, dy_est = %f, error = %f\n',...
%                 iter,position_est,dx_est(position_est),dy_est(position_est),err);
        end
        iter = iter + 1;

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
    cord_est = cord_fixed_est + [dx_est(position_est);dy_est(position_est)];
    h_est_N = h_est_mat(position_est, :).';
    
    supp = position_est : P.Q : position_est+(P.N-1)*P.Q;
    Atmp = cell2mat(AOAs);
    AOAs_est = Atmp(supp);
%     fprintf('\ncord_fixed_est = [%f, %f]\n',cord_fixed_est(1), cord_fixed_est(2));
end

