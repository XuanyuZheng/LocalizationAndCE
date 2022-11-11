function [cord_est] = proposed_estPHN_DynGrid(y,SNRs,P,PHN_deg)
%proposed_estPHN_DynGrid 
%Estimate PHN, Constrained formulation, Dynamic Grid Update

    % Gamma dist parameters
    beta = 0.0001;
    alpha = 0.01;
    tx_powers = (10.^(SNRs/10));
%     sigma_x = sqrt(var_x);
    % get y for each BS
    y_cell = zeros(P.Nr,P.N);
    for n = 1:P.N
        y_cell(:,n) = y((n-1)*P.Nr + 1 : n*P.Nr);
    end

    % Initialization of PHN_est, h_est and dx_est and dy_est
    plus_minus = unidrnd(2,P.N*P.Q,1);
    for i = 1:P.N*P.Q
       if  plus_minus(i) == 2
           plus_minus(i) = -1;
       end
    end
    h_est = plus_minus*(1/P.Q); % initialize to mean of h_tilde
    % h_est = ones(N*Q,1);      % % initialize to zeros
    PHN_est = zeros(P.Nr,P.N);   % or just rand
    dx_est = zeros(P.Q,1);
    dy_est = zeros(P.Q,1);
    % Finished Initialization


    %Algorithm specifications
    thres_fix = 1e-6;               % stopping criteria for h for fixed grid
    thres = 1e-4;                   % stopping criteria for h after grid refinement
    err = Inf;                      % First error
    max_iter = 30;                  % max number of iteration
    bound = (  1 + 2 * ((10.^(max(SNRs)/10)) * (1-cos(deg2rad(1*PHN_deg))))  );
    % bound = norm(y-Omega_fix*h_true)^2/P.N/P.Nr
    % bound = 1
    fix_converged = 0;
    iter = 1;
    % update sensing matrix Omega & update steerring matrix Ar
    [Omega,~,~] = ...
    getSpsRecry(dx_est,dy_est,SNRs,PHN_est,P);
    obj_value = norm(y-Omega*h_est)^2;
%     obj_value = obj_lik(y, SNRs, PHN_est, dx_est, dy_est, h_est, P);
    while iter <= max_iter && err >= thres

%         obj_value_old = obj_value;
        % Fix Lambda, dx, dy, solve h
        h_est_old = h_est;

        if iter == 1 %&& SNR <= 10
            inner_it = 1; %% change to 5?
        else
            inner_it = 1;
        end
        for h_iter = 1:inner_it
            
            % update weights
            lambda = zeros(P.Q,1);
            for q = 1:P.Q
                h_bar_sum_over_n = 0;
                for n = 1:P.N
                    h_bar_sum_over_n = h_bar_sum_over_n + norm(h_est(   (n-1)*P.Q+q   ),2)^2;       %Here use prior norm????????????
                end
                lambda(q) = (alpha+P.N)/(beta+sqrt(h_bar_sum_over_n));
            end
            %solve h
            % do variable substution, s.t. x(n,q) = lambda(q)*h(n,q)
            A = zeros(P.N*P.Nr,P.N*P.Q);
            for n1 = 1:P.N
               for q1 = 1:P.Q
                  A(:,(n1-1)*P.Q+q1) = Omega(:,(n1-1)*P.Q+q1)/(lambda(q1));
               end
            end
            y1 = y;
            cvx_begin quiet
                cvx_solver mosek
                variable x1(P.N,P.Q) complex
    %             minimize( norm(x1,2)  )
                minimize( sum(norms(x1))  )
                    subject to
                    sum_square_abs(y1-A*vec(transpose(x1))) <= (1.01*P.N*P.Nr*bound); %#ok<VUNUS>
            cvx_end
            %toc
            for n1 = 1:P.N
               for q1 = 1:P.Q
                  h_est((n1-1)*P.Q+q1) = x1(n1,q1)/(lambda(q1));
               end
            end
        end
        
%         bound = bound/2;
%         if bound < 1.1
%             bound = 1.1;
%         end
        [~,Ar,~] = ...
        getSpsRecry(dx_est,dy_est,SNRs,PHN_est,P);
%         if isnan(h_est(1))
%             h_est(1)
%             break;
%         end
        % Finished updating h_est

        % Fix h, dx, dy, solve Lambda
        h_est_cell = cell(P.N,1);
        for i = 1:P.N
            SNR = SNRs(i);
            tx_power = (10^(SNR/10));
            h_tilde_n = h_est((i-1)*P.Q+1 : i*P.Q);
            h_est_cell{i} = h_tilde_n;
            y_tilde_n = P.x*sqrt(tx_power)*Ar{i}*h_tilde_n;
           for j = 1:P.Nr
%                 PHN_est(j,i) = angle(y((i-1)*P.Nr+j)) - angle(y_tilde_n(j));
                PHN_est(j,i) = 0;
           end
        end

        % update sensing matrix Omega & update steerring matrix Ar
        [Omega,Ar,AOA] = ...
        getSpsRecry(dx_est,dy_est,SNRs,PHN_est,P);

        if fix_converged == 0
            err_fixed = norm(h_est-h_est_old)/norm(h_est_old);
            if err_fixed < thres_fix
                fix_converged = 1;
            end
        else
            err = norm(h_est-h_est_old)/norm(h_est_old);
        end
       


        % Fix Lambda, h solve dx_est, dy_est
        if fix_converged == 1
%             obj_value_old = obj_value;

            dx_grad = zeros(P.Q,1);
            dy_grad = zeros(P.Q,1);
            
            P_nom = tx_powers;
%             P_nom = var_x/P.var_n;      %normalized tx power

            % Construct [yn]-q
%             Y_sub = zeros(P.Nr,P.N,P.Q);
%             for n = 1:P.N
%                 Xn = sqrt(P_nom)*P.x*rhos(n)*diag(exp(1i*PHN_est(:,n)));
%                 for q = 1:P.Q
%                     Y_sub(:,n,q) = y_cell(:,n) - ( Xn*Ar{n}*h_est_cell{n} -  Xn*Ar{n}(:,q)*h_est((n-1)*P.Q+q));
%                 end
%             end
            
            for i = 1:P.Q
                dx_grad_this = 0;
%                 if i == 234
%                     keyboard;
%                 end
                for n = 1:P.N

%                     dtheta_ni_dx_i = (P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))/(((P.qs{i}(1)+dx_est(i))-P.BS_Pos(1,n))^2*...
%                         (1+((P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))/(P.BS_Pos(1,n)-(P.qs{i}(1)+dx_est(i))))^2));
% 
%                     a_n_i = Ar{n}(:,i);
%                     % a_n_i_tilde = a_n_i.*(Ant_pos{n}*[-sin(AOA{n}(i));cos(AOA{n}(i))]);
% 
%                     yn_i = Y_sub(:,n,i);
%                     h_n_i = h_est((n-1)*P.Q+i);
% 
%                     a_nq_diff_x = dtheta_ni_dx_i*2*pi*1i/P.lambda_c*a_n_i.*(P.antenna_pos'*[-sin(AOA{n}(i));cos(AOA{n}(i))]);
% 
%                     dx_grad_this = dx_grad_this + 2*real(a_nq_diff_x'*(Xn'*Xn)*a_n_i)*abs(h_n_i)^2 ...
%                         -2*real(a_nq_diff_x'*Xn'*h_n_i'*yn_i);
% 
%                 end
%                 dx_grad(i) = dx_grad_this;
                    dtheta_ni_dx_i = (P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))/(((P.qs{i}(1)+dx_est(i))-P.BS_Pos(1,n))^2*...
                        (1+((P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))/(P.BS_Pos(1,n)-(P.qs{i}(1)+dx_est(i))))^2));
                    a_n_i = Ar{n}(:,i);
                    a_n_i_tilde = a_n_i.*(P.antenna_pos'*[-sin(AOA{n}(i));cos(AOA{n}(i))]);

                    Lmabda_n = diag(exp(1i*PHN_est(:,n)));
                    h_n_i = h_est((n-1)*P.Q+i);
        %             if n == 4 && i == 46
        %                 h_n_i = h_n_i';
        %             end

                    dx_grad_this = dx_grad_this + 4*pi/P.lambda_c*dtheta_ni_dx_i*(...
                        -sqrt(P_nom(n))*P.x*real(1i*h_n_i*(y_cell(:,n))' * Lmabda_n * a_n_i_tilde) ...
                        + P_nom(n)*abs(P.x)^2*real(1i*h_n_i*h_est_cell{n}'*Ar{n}'*a_n_i_tilde)...
                        );
                end
                dx_grad(i) = dx_grad_this;
            end
            % back tracking
            step_size = 1;
            alpha_bt = 0.1;
            beta_bt = 0.3;
            
            obj_this_x = obj_lik(y,SNRs,PHN_est,dx_est,dy_est,h_est,P);
            while 1
                obj_dyn_x = obj_lik(y,SNRs,PHN_est,dx_est+step_size*(-dx_grad),dy_est,h_est, P);
                if obj_dyn_x < obj_this_x + alpha_bt*step_size*dx_grad'*(-dx_grad)
                    break
                else
                    step_size = beta_bt*step_size;
                end
            end
            dx_est = dx_est + step_size * (-dx_grad);
            
            
            % update sensing matrix Omega & update steerring matrix Ar
            [~,Ar,AOA] = ...
            getSpsRecry(dx_est,dy_est,SNRs,PHN_est,P);
        
            for i = 1:P.Q
                dy_grad_this = 0;
%                 if i == 157
%                     keyboard;
%                 end
                for n = 1:P.N
%                     dtheta_ni_dy_i = 1/(((P.qs{i}(1)+dx_est(i))-P.BS_Pos(1,n))*(1+((P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))...
%                         /(P.BS_Pos(1,n)-(P.qs{i}(1)+dx_est(i))))^2));
%                     a_n_i = Ar{n}(:,i);
%                     % a_n_i_tilde = a_n_i.*(Ant_pos{n}*[-sin(AOA{n}(i));cos(AOA{n}(i))]);
% 
%                     yn_i = Y_sub(:,n,i);
%                     h_n_i = h_est((n-1)*P.Q+i);
% 
%                     a_nq_diff_y = dtheta_ni_dy_i*2*pi*1i/P.lambda_c*a_n_i.*(P.antenna_pos'*[-sin(AOA{n}(i));cos(AOA{n}(i))]);
% 
%                     dy_grad_this = dy_grad_this + 2*real(a_nq_diff_y'*(Xn'*Xn)*a_n_i)*abs(h_n_i)^2 ...
%                         -2*real(a_nq_diff_y'*Xn'*h_n_i'*yn_i);
%                 end
%                 dy_grad(i) = dy_grad_this;
                    dtheta_ni_dy_i = 1/(((P.qs{i}(1)+dx_est(i))-P.BS_Pos(1,n))*(1+((P.BS_Pos(2,n)-(P.qs{i}(2)+dy_est(i)))...
                        /(P.BS_Pos(1,n)-(P.qs{i}(1)+dx_est(i))))^2));
                    a_n_i = Ar{n}(:,i);
                    a_n_i_tilde = a_n_i.*(P.antenna_pos'*[-sin(AOA{n}(i));cos(AOA{n}(i))]);
                    Lmabda_n = diag(exp(1i*PHN_est(:,n)));
                    h_n_i = h_est((n-1)*P.Q+i);
                    dy_grad_this = dy_grad_this + 4*pi/P.lambda_c*dtheta_ni_dy_i*(...
                        -sqrt(P_nom(n))*P.x*real(1i*h_n_i*(y_cell(:,n))' * Lmabda_n * a_n_i_tilde) ...
                        + P_nom(n)*abs(P.x)^2*real(1i*h_n_i*h_est_cell{n}'*Ar{n}'*a_n_i_tilde)...
                        );   
                end
                dy_grad(i) = dy_grad_this;
            end
            % Finished computing gradient
            
%             % fixed step size update       
%             dx_grad = dx_grad/norm(dx_grad);
%             dy_grad = dy_grad/norm(dy_grad);


%             dx_est = dx_est - dx_grad * 0.1;
%             dy_est = dy_est - dy_grad * 0.1; 
  
            step_size = 1;
            obj_this_y = obj_lik(y,SNRs,PHN_est,dx_est,dy_est,h_est,P);
            while 1
                obj_dyn_y = obj_lik(y,SNRs,PHN_est,dx_est,dy_est + step_size*(-dy_grad),h_est, P);
                if obj_dyn_y < obj_this_y + alpha_bt*step_size*dy_grad'*(-dy_grad)
                    break
                else
                    step_size = beta_bt*step_size;
                end
            end
            dy_est = dy_est + step_size * (-dy_grad);
            
            
            % update sensing matrix Omega & update steerring matrix Ar
            [Omega,~,~] = ...
            getSpsRecry(dx_est,dy_est,SNRs,PHN_est,P);
        
%             obj_value = norm(y-Omega*h_est)^2;
            obj_value = obj_lik(y,SNRs,PHN_est,dx_est,dy_est,h_est,P);
%             vpa(obj_value)
%             if (obj_value > obj_value_old)
%                 break;
%             end

            % update sensing matrix Omega & update steerring matrix Ar
%             [Omega,~,~] = ...
%             getSpsRecry(dx_est,dy_est,sigma_x,rhos,PHN_est,P);
        

%             err = norm(h_est-h_est_old)/norm(h_est_old);
        end
        vpa(obj_value)
%         if obj_value > obj_value_old
%             break;
%         end
        
    
        % other analysis
        
        
%         % trascking PHN error
%         PHN_error = mod(mod(PHN,2*pi)-mod(PHN_est,2*pi),2*pi);
%         for i = 1:P.Nr
%            for j = 1:P.N
%                if(PHN_error(i,j))>pi
%                    PHN_error(i,j) = PHN_error(i,j)-2*pi;
%                end
%                if(PHN_error(i,j))<-pi
%                    PHN_error(i,j) = PHN_error(i,j)+2*pi;
%                end
%            end
%         end
%         PHN_MSE = zeros(P.N,1);                                               
%         for i = 1:P.N
%             PHN_MSE(i) = norm(PHN_error(:,i),2)^2/P.Nr;
%         end
%         PHN_MSE_all = sum(PHN_MSE)/P.N;                                   % PHN mean square error, recorded!

%         % objective function

%         h_est_mat_tr = transpose(h_est_mat);
%         h_q_abs = transpose(sum(abs(h_est_mat_tr)));
%         log_ph_sum = 0;     %log(p(h))
%         for q = 1:P.Q
%             log_ph_sum = log_ph_sum + ...
%                 -(P.N+alpha)*log((1+h_q_abs(q)));
%         end
%         sigma = 1;
%         f = sigma^-1*(norm(y-Omega*h_est))^2-log_ph_sum + 1/(2*sigma^2)*norm(PHN_est,2)^2;                %should decrease
%         log_ph_sum_tr = -(P.N+alpha)*log((1+norm(h_tilde,1)));
%         f_tr = sigma^-1*(norm(y-Omega_tr*h_true))^2-log_ph_sum_tr + 1/(2*sigma^2)*norm(PHN,2)^2;
%         f_err = norm(f-f_tr)/norm(f_tr);

    %     percentage = ((m-1)*var_x_dBms_num*test_num+(k-1)*test_num+test_iter)/(PHN_num*var_x_dBms_num*test_num);
%         fprintf('PHN_deg=%d, TxPower=%d, iter=%d, obj=%f, obj_tr=%f, obj_err=%f, err=%f, PHN_err=%f\n',...
%             PHN_deg,var_x_dBm,iter,f,f_tr,f_err,err,PHN_MSE_all);
        h_est_mat = zeros(P.Q,P.N);
        for i = 1:P.N
           for j = 1:P.Q
               h_est_mat(j,i) = h_est((i-1)*P.Q+j);
           end
        end
        h_est_mat_sum = sum(abs(h_est_mat'));
        plot((h_est_mat_sum'))
        position_est = find(h_est_mat_sum==max(h_est_mat_sum));          % estimated position, recorded!
        fprintf('iter=%d, position_est = %d, dx_est = %f, dy_est = %f, obj_value = %f, error = %f, error_fixed = %f\n',...
            iter,position_est,dx_est(position_est),dy_est(position_est),obj_value,err,err_fixed);
        iter = iter + 1;
    end

    h_est_mat_sum = sum(abs(h_est_mat'));
    plot((h_est_mat_sum'))
    position_est = find(h_est_mat_sum==max(h_est_mat_sum));          % estimated position, recorded!

    cord_fixed_est = [(mod(position_est,P.grid_side)-1)*P.grid_length+P.grid_length/2-P.grid_side*P.grid_length/2;...
    -((ceil(position_est/P.grid_side)-1)*P.grid_length+P.grid_length/2)+P.grid_side*P.grid_length/2];

    if mod(position_est,P.grid_side) == 0         %deal with mod = 0 but 
        cord_fixed_est(1) = (P.grid_side-1)*P.grid_length+P.grid_length/2-P.grid_side*P.grid_length/2;
    end
    fprintf('\ncord_fixed_est = [%f, %f],           alpha = %f\n',cord_fixed_est(1), cord_fixed_est(2), alpha);
    cord_est = cord_fixed_est + [dx_est(position_est);dy_est(position_est)];
end

