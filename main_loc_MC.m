%% Let's get started
clear;
close;
tic;
% parpool('local',48);               % Use parallel computting tool
%%
P = setUpSystem();
PHN_degs = 5;
PHN_num = length(PHN_degs);
var_x_dBms = 30;
var_x_dBms_num = length(var_x_dBms);

max_test = 1;


err_loc_test_PHN_Vdb = zeros(max_test,PHN_num,var_x_dBms_num);
err_ch_amp_test_PHN_Vdb = zeros(P.N,max_test,PHN_num,var_x_dBms_num);
err_ch_test_PHN_Vdb = zeros(P.N,max_test,PHN_num,var_x_dBms_num);
err_AOA_test_PHN_Vdb = zeros(P.N,max_test,PHN_num,var_x_dBms_num);

err_loc_test_PHN_Vdb_unc = zeros(max_test,PHN_num,var_x_dBms_num);
err_ch_amp_test_PHN_Vdb_unc = zeros(P.N,max_test,PHN_num,var_x_dBms_num);
err_ch_test_PHN_Vdb_unc = zeros(P.N,max_test,PHN_num,var_x_dBms_num);
err_AOA_test_PHN_Vdb_unc = zeros(P.N,max_test,PHN_num,var_x_dBms_num);

err_loc_test_PHN_Vdb_NoPHN = zeros(max_test,PHN_num,var_x_dBms_num);
err_ch_amp_test_PHN_Vdb_NoPHN = zeros(P.N,max_test,PHN_num,var_x_dBms_num);
err_ch_test_PHN_Vdb_NoPHN = zeros(P.N,max_test,PHN_num,var_x_dBms_num);
err_AOA_test_PHN_Vdb_NoPHN = zeros(P.N,max_test,PHN_num,var_x_dBms_num);

Wrong_num_PHN_Vdb = zeros(PHN_num, var_x_dBms_num);
small_h_num_PHN_Vdb = zeros(PHN_num, var_x_dBms_num);
small_h_num_PHN_Vdb_NoPHN = zeros(PHN_num, var_x_dBms_num);


for v = 1:var_x_dBms_num
    
    var_x_dBm = var_x_dBms(v);



    err_loc_test_PHN = zeros(max_test,PHN_num);
    err_ch_amp_test_PHN = zeros(P.N,max_test,PHN_num);
    err_ch_test_PHN = zeros(P.N,max_test,PHN_num);
    err_AOA_test_PHN = zeros(P.N,max_test,PHN_num);
    
    err_loc_test_PHN_unc = zeros(max_test,PHN_num);
    err_ch_amp_test_PHN_unc = zeros(P.N,max_test,PHN_num);
    err_ch_test_PHN_unc = zeros(P.N,max_test,PHN_num);
    err_AOA_test_PHN_unc = zeros(P.N,max_test,PHN_num);


    err_loc_test_PHN_NoPHN = zeros(max_test,PHN_num);
    err_ch_amp_test_PHN_NoPHN = zeros(P.N,max_test,PHN_num);
    err_ch_test_PHN_NoPHN = zeros(P.N,max_test,PHN_num);
    err_AOA_test_PHN_NoPHN = zeros(P.N,max_test,PHN_num);
    
    Wrong_num_PHN = zeros(PHN_num,1);
    small_h_num_PHN = zeros(PHN_num,1);
    small_h_num_PHN_NoPHN = zeros(PHN_num,1);

    for s = 1:PHN_num
        PHN_deg = PHN_degs(s);


        err_loc_test = zeros(max_test,1);
        err_ch_amp_test = zeros(P.N,max_test);
        err_ch_test = zeros(P.N,max_test);
%         err_PHN_test = zeros(P.N*P.Nr,max_test);
%         mse_PHN_test = zeros(max_test,1);
        err_AOA_test = zeros(P.N,max_test);
%         err_PHN_re_test = zeros(P.N*P.Nr,max_test);
%         mse_PHN_re_test = zeros(max_test,1);

        err_loc_test_unc = zeros(max_test,1);
        err_ch_amp_test_unc = zeros(P.N,max_test);
        err_ch_test_unc = zeros(P.N,max_test);
%         err_PHN_test = zeros(P.N*P.Nr,max_test);
%         mse_PHN_test = zeros(max_test,1);
        err_AOA_test_unc = zeros(P.N,max_test);

        err_loc_test_NoPHN = zeros(max_test,1);
        err_ch_amp_test_NoPHN = zeros(P.N,max_test);
        err_ch_test_NoPHN = zeros(P.N,max_test);
%         err_PHN_test_NoPHN = zeros(P.N*P.Nr,max_test);
%         mse_PHN_test_NoPHN = zeros(max_test,1);
        err_AOA_test_NoPHN = zeros(P.N,max_test);
        
        Wrong_num = 0;
        small_h_num = 0;
        small_h_num_No = 0;

%         test = 1;
        for test = 1: max_test


            %% System parameters


            % position of the user
            row_index = randi([1 P.grid_side]);                         % pick one grid from x axis
            col_index = randi([1 P.grid_side]);                         % pick one grid from y axis
            MU_Pos_fixed = P.qs_mesh{row_index,col_index};              % pick out the rand grid
            qs_mat = [P.qs{:}];                                         % 2 x 400 matrix containing the grid coordinate (left to right, up to down)
            pos_index_tr = (row_index-1)*P.grid_side + col_index;       % real index of the grid where the user is in
            x_off_tr = unifrnd(-P.grid_length/2,P.grid_length/2);       % x axis off-grid
            y_off_tr = unifrnd(-P.grid_length/2,P.grid_length/2);       % y axis off-grid
            % x_off_tr = 0;      % x axis off-grid
            % y_off_tr = 0;      % y axis off-grid
            MU_pos_tr = MU_Pos_fixed + [x_off_tr;y_off_tr];             % true coordinate of the user

            % position of the scatterers
            scat_pos_tr = [[10;175] [-20;-225]];
            % scat_pos_tr = [[-200;200] [200;200] [-200;-200] [200;-200]];

            effect = [1 1 0 0;0 0 1 1].';
%             effect = [1 1 1 1;1 1 1 1].';
            % effect = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1].';

            K = size(scat_pos_tr, 2);

            %% compute pathloss/receive sinr
            pl_exp_LoS = 2;     % LoS pathloss exponent
            pl_exp_NLoS = 2.8;    % NLoS pathloss exponent
            rhos_LoS = getPathLoss(MU_pos_tr, P, pl_exp_LoS);   % LoS pathloss
            rhos_NLoS = zeros(P.N, K);                   % NLoS pathloss
            for k = 1:K
                rhos_NLoS(:,k) = getPathLossNLoS(MU_pos_tr,scat_pos_tr(:,k), P, pl_exp_NLoS);
            end

            rhos_NLoS_sig = rhos_NLoS.*effect;  %transmitted path loss

            var_x = 10^(var_x_dBm/10)/1000;                     % tx power/variance in w
            SNRs_LoS = 10*log10(var_x.*rhos_LoS.^2./P.var_n);           % LoS receive SNRs in N BSs
            SNRs_NLoS = 10*log10(var_x.*rhos_NLoS.^2./P.var_n);           % LoS receive SNRs in N BSs
            tx_power_LoS = (10.^(SNRs_LoS/10));
            tx_power_NLoS = (10.^(SNRs_NLoS/10));

            SNRs_NLoS_sig = 10*log10(var_x.*rhos_NLoS_sig.^2./P.var_n);           % LoS receive SNRs in N BSs
            tx_power_NLoS_sig = (10.^(SNRs_NLoS_sig/10));

            % PHN_tr is the true PHN, w_N_tr is the true LoS channel gains, v_NK_tr is the true NLoS channel gains
            % h_tr is the true overall channel, CovMat is the covariance matrix before multiplying PHN_deg
            [y,PHN_tr,w_N_tr,v_NK_tr,h_tr,AoA_NLoS_tr,CovMat] = ...
                getAllReceivedSigs(MU_pos_tr,scat_pos_tr,SNRs_LoS,SNRs_NLoS_sig,PHN_deg,P);
            PHN_CovMat_n = deg2rad(PHN_deg)^2*CovMat;
            PHN_CovMat = [];
            for n = 1:P.N
                PHN_CovMat = blkdiag(PHN_CovMat,PHN_CovMat_n);
            end
            %% groud truth
            PHN_tr;     % true phn
            w_N_tr;     % true LoS channel gains, N x 1
            v_NK_tr;    % true NLoS channel gains, N x K
            h_tr;       % true overall channels, Nr x N
            AoA_NLoS_tr;% true AOAs of NLoS paths, N x K (-90~270)
            w_pl_N_tr = sqrt(tx_power_LoS).*w_N_tr;
            v_pl_NK_tr = sqrt(tx_power_NLoS).*v_NK_tr;



            pos_index_tr;                       % true support for LoS
            pos_index_NLoS_tr = zeros(P.N,K);   % true support for NLoS (0~M for each NLoS path)
            for n=1:P.N
                for k = 1:K
                    for ind = 1:P.M - 1
                        pos_index_NLoS_tr(n,k) = P.M;
                        if AoA_NLoS_tr(n,k) >= P.phi(ind) ...
                                && AoA_NLoS_tr(n,k) < P.phi(ind+1)
                            pos_index_NLoS_tr(n,k) = ind;
                            break;
                        end
                    end
                end
            end
            % pos_index_NLoS_tr = pos_index_NLoS_tr.*effect;
            % actual sparse recovery
            dx = zeros(P.Q,1);
            dy = zeros(P.Q,1);
            dx(pos_index_tr) = x_off_tr;
            dy(pos_index_tr) = y_off_tr;

            dphi = zeros(P.M,P.N);
            for n = 1:P.N
                for k = 1:K
                    this_supp = pos_index_NLoS_tr(n,k);
                    if effect(n,k) ~= 0
                        dphi(this_supp,n) = ...
                            AoA_NLoS_tr(n,k) - P.phi_mesh(this_supp);
                    end
                end
            end

            % vectorized true channel
            w_true = zeros(P.N*P.Q,1);
            supp_w = pos_index_tr : P.Q : pos_index_tr+(P.N-1)*P.Q;
            w_true(supp_w) = w_N_tr;

            v_true = zeros(P.N*P.M,1);
            v_pl_true = zeros(P.N*P.M,1);
            supp_v = zeros(P.N*K,1);
            for n = 1:P.N
                supp_v((n-1)*K + 1:n*K ) = ((pos_index_NLoS_tr(n,:).')+(n-1)*P.M);
            end

            supp_v_effect = zeros(P.N*K,1);
            for n = 1:P.N
                supp_v_effect((n-1)*K + 1:n*K ) = ((pos_index_NLoS_tr(n,:).')+(n-1)*P.M).*(effect(n,:).');
            end
            supp_v_effect = supp_v_effect(supp_v_effect~=0);

            v_NK_tr_effect = v_NK_tr.*effect;
            v_true(supp_v) = vec((v_NK_tr_effect).');

            v_pl_NK_tr_effect = v_pl_NK_tr.*effect;
            v_pl_true(supp_v) = vec(v_pl_NK_tr_effect.');

            % get true sparse recovery
            [Omega_tr,Ars_tr,AOAs_qs_tr] = ...
                getSpsRecry_LoS(dx,dy,SNRs_LoS,P);  % involved pathloss in Phi

            [Phi_tr,Brs_tr,AOAs_phis_tr] = ...
                getSpsRecry_NLoS(dphi,SNRs_NLoS,P);           % didn't involve pathloss in Phi
%             y_err = y-Omega_tr*w_true -Phi_tr*v_true;

        %% Start the algorithm of localization with PHN with NLoS and dynamic grid update
%         [cord_est,h_est,h_est_source] = loc_SRED_SPI_dyn(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P);
% 
%         h_err = norm(h_est - h_tr)/norm(h_tr);
%         h_err_source = norm(h_est_source - h_tr)/norm(h_tr);
%         loc_err = norm(cord_est - MU_pos_tr)
%         20*log10(h_err)
%         20*log10(h_err_source)

        %% ignore PHN, uncinstrained
%         [cord_est_NoPHN,h_est_NoPHN,h_est_source_NoPHN] = loc_SRED_SPI_dyn_NoPHN(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P);
% 
%         h_err_NoPHN = norm(h_est_NoPHN - h_tr)/norm(h_tr);
%         h_err_source_NoPHN = norm(h_est_source_NoPHN - h_tr)/norm(h_tr);
%         loc_err_NoPHN = norm(cord_est_NoPHN - MU_pos_tr)
%         20*log10(h_err_NoPHN)
%         20*log10(h_err_source_NoPHN)

            %% Start of algorithm

            try
                [cord_est,h_est,h_est_source] = loc_SRED_SPI_dyn(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P);
                [cord_est_unc,h_est_unc,h_est_source_unc] = loc_SRED_SPI_dyn_unc(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P);
                [cord_est_NoPHN,h_est_NoPHN,h_est_source_NoPHN] = loc_SRED_SPI_dyn_NoPHN(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P);
%                 [cord_est, h_est_N, AOAs_est] = proposed_new_likelihood_Dyn(y,SNRs,PHN_deg,PHN_CovMat,P);
%                 [cord_est_unc, h_est_N_unc, AOAs_est_unc] = proposed_new_likelihood_uncons_Dyn(y,SNRs,PHN_deg,PHN_CovMat,P);
%                 [cord_est_NoPHN, h_est_N_NoPHN, AOAs_est_NoPHN] = proposed_new_likelihood_NoPHN_Dyn(y,SNRs,PHN_deg,PHN_CovMat,P);
                
            catch
    %             test = test - 1;
                fprintf('Wrong')
                Wrong_num = Wrong_num + 1;
                continue
            end
    %         test

            
            err_loc = norm(cord_est-MU_pos_tr);	err_loc_unc = norm(cord_est_unc-MU_pos_tr);	err_loc_NoPHN = norm(cord_est_NoPHN-MU_pos_tr);% 1 x 1
%             err_ch_amp = abs((abs(h_est_N) - abs(h))) ./ abs(h);    err_ch_amp_unc = abs((abs(h_est_N_unc) - abs(h))) ./ abs(h);	err_ch_amp_NoPHN = abs((abs(h_est_N_NoPHN) - abs(h))) ./ abs(h);    % 4 x 1
%             err_ch = abs(h_est_N - h) ./ abs(h);    err_ch_unc = abs(h_est_N_unc - h) ./ abs(h);	err_ch_NoPHN = abs(h_est_N_NoPHN - h) ./ abs(h);                    % 4 x 1
    %         err_PHN = wrapToPi(PHN(:)-PHN_est(:));                  err_PHN_NoPHN = wrapToPi(PHN(:)-PHN_est_NoPHN(:));                  % 80 x 1
    %         mse_PHN = sqrt(sum(err_PHN.^2)/(P.N*P.Nr));             mse_PHN_NoPHN = sqrt(sum(err_PHN_NoPHN.^2)/(P.N*P.Nr));             % 1 x 1
%             err_AOA = wrapToPi(AOAs_est-AOAs_tr);   err_AOA_unc = wrapToPi(AOAs_est_unc-AOAs_tr);	err_AOA_NoPHN = wrapToPi(AOAs_est_NoPHN-AOAs_tr);
    %         
    %         h_ph = zeros(P.Nr,P.N);
    %         h_ph_est = zeros(P.Nr,P.N);
    %         for i = 1:P.N
    %             h_ph(:,i) = ones(P.Nr,1) * angle(h(i));
    %             h_ph_est(:,i) = ones(P.Nr,1) * angle(h_est_N(i));
    %         end
    %         h_ph_vec = vec(h_ph);
    %         h_ph_diag = diag(h_ph_vec);
    %         PHN_re = h_ph_diag * PHN(:);
    % 
    %         h_ph_vec_est = vec(h_ph_est);
    %         h_ph_diag_est = diag(h_ph_vec_est);
    %         PHN_re_est = h_ph_diag_est * PHN_est(:);
    % 
    %         err_PHN_re = wrapToPi(PHN_re(:)-PHN_re_est(:));               % 80 x 1
    %         mse_PHN_re = sqrt(sum(err_PHN_re.^2)/(P.N*P.Nr));             % 1 x 1
    % 
    % 
    % 
            err_loc_test(test) = err_loc;
%             err_ch_amp_test(:,test) = err_ch_amp;
%             err_ch_test(:,test) = err_ch;
    %         err_PHN_test(:,test) = err_PHN;
    %         mse_PHN_test(test) = mse_PHN;
%             err_AOA_test(:,test) = err_AOA;
    %         err_PHN_re_test(:,test) = err_PHN_re;
    %         mse_PHN_re_test(test) = mse_PHN_re;
    
            err_loc_test_unc(test) = err_loc_unc;
%             err_ch_amp_test_unc(:,test) = err_ch_amp_unc;
%             err_ch_test_unc(:,test) = err_ch_unc;
    %         err_PHN_test(:,test) = err_PHN;
    %         mse_PHN_test(test) = mse_PHN;
%             err_AOA_test_unc(:,test) = err_AOA_unc;
    %         
            err_loc_test_NoPHN(test) = err_loc_NoPHN;
%             err_ch_amp_test_NoPHN(:,test) = err_ch_amp_NoPHN;
%             err_ch_test_NoPHN(:,test) = err_ch_NoPHN;
%     %         err_PHN_test_NoPHN(:,test) = err_PHN_NoPHN;
%     %         mse_PHN_test_NoPHN(test) = mse_PHN_NoPHN;
%             err_AOA_test_NoPHN(:,test) = err_AOA_NoPHN;
%     %         
%             err_AOA_this = rad2deg(sqrt(sum(err_AOA.^2)/P.N));
%             err_AOA_this_NoPHN = rad2deg(sqrt(sum(err_AOA_NoPHN.^2)/P.N));
            fprintf('Vdb = %d, PHN_deg = %d, test = %d, err_loc = %f, err_loc_unc = %f, err_loc_No = %f\n',...
                var_x_dBm, PHN_deg, test, err_loc, err_loc_unc, err_loc_NoPHN);


%             test = test + 1;

        end
        Wrong_num_PHN(s) = Wrong_num;
        small_h_num_PHN(s) = small_h_num;
        small_h_num_PHN_NoPHN(s) = small_h_num_No;
        
        err_loc_test_PHN(:,s) = err_loc_test;
        err_ch_amp_test_PHN(:,:,s) = err_ch_amp_test;
        err_ch_test_PHN(:,:,s) = err_ch_test;
    %     err_PHN_test_PHN(:,:,s) = err_PHN_test;
    %     mse_PHN_test_PHN(:,s) = mse_PHN_test;
        err_AOA_test_PHN(:,:,s) = err_AOA_test;
    %     err_PHN_re_test_PHN(:,:,s) = err_PHN_re_test;
    %     mse_PHN_re_test_PHN(:,s) = mse_PHN_re_test;
    
        err_loc_test_PHN_unc(:,s) = err_loc_test_unc;
        err_ch_amp_test_PHN_unc(:,:,s) = err_ch_amp_test_unc;
        err_ch_test_PHN_unc(:,:,s) = err_ch_test_unc;
    %     err_PHN_test_PHN(:,:,s) = err_PHN_test;
    %     mse_PHN_test_PHN(:,s) = mse_PHN_test;
        err_AOA_test_PHN_unc(:,:,s) = err_AOA_test_unc;
    %     
        err_loc_test_PHN_NoPHN(:,s) = err_loc_test_NoPHN;
        err_ch_amp_test_PHN_NoPHN(:,:,s) = err_ch_amp_test_NoPHN;
        err_ch_test_PHN_NoPHN(:,:,s) = err_ch_test_NoPHN;
%     %     err_PHN_test_PHN_NoPHN(:,:,s) = err_PHN_test_NoPHN;
%     %     mse_PHN_test_PHN_NoPHN(:,s) = mse_PHN_test_NoPHN;
        err_AOA_test_PHN_NoPHN(:,:,s) = err_AOA_test_NoPHN;

    end
    
    Wrong_num_PHN_Vdb(:,v) = Wrong_num_PHN;
    small_h_num_PHN_Vdb(:,v) = small_h_num_PHN;
    small_h_num_PHN_Vdb_NoPHN(:,v) = small_h_num_PHN_NoPHN;
    
    err_loc_test_PHN_Vdb(:,:,v) = err_loc_test_PHN;
    err_ch_amp_test_PHN_Vdb(:,:,:,v) = err_ch_amp_test_PHN;
    err_ch_test_PHN_Vdb(:,:,:,v) = err_ch_test_PHN;
    err_AOA_test_PHN_Vdb(:,:,:,v) = err_AOA_test_PHN;
    
    err_loc_test_PHN_Vdb_unc(:,:,v) = err_loc_test_PHN_unc;
    err_ch_amp_test_PHN_Vdb_unc(:,:,:,v) = err_ch_amp_test_PHN_unc;
    err_ch_test_PHN_Vdb_unc(:,:,:,v) = err_ch_test_PHN_unc;
    err_AOA_test_PHN_Vdb_unc(:,:,:,v) = err_AOA_test_PHN_unc;

    err_loc_test_PHN_Vdb_NoPHN(:,:,v) = err_loc_test_PHN_NoPHN;
    err_ch_amp_test_PHN_Vdb_NoPHN(:,:,:,v) = err_ch_amp_test_PHN_NoPHN;
    err_ch_test_PHN_Vdb_NoPHN(:,:,:,v) = err_ch_test_PHN_NoPHN;
    err_AOA_test_PHN_Vdb_NoPHN(:,:,:,v) = err_AOA_test_PHN_NoPHN;
    
    save('err_loc_test_PHN_Vdb.mat','err_loc_test_PHN_Vdb');
    save('err_ch_amp_test_PHN_Vdb.mat','err_ch_amp_test_PHN_Vdb');
    save('err_ch_test_PHN_Vdb.mat','err_ch_test_PHN_Vdb');
    save('err_AOA_test_PHN_Vdb.mat','err_AOA_test_PHN_Vdb');
    
    save('err_loc_test_PHN_Vdb_unc.mat','err_loc_test_PHN_Vdb_unc');
    save('err_ch_amp_test_PHN_Vdb_unc.mat','err_ch_amp_test_PHN_Vdb_unc');
    save('err_ch_test_PHN_Vdb_unc.mat','err_ch_test_PHN_Vdb_unc');
    save('err_AOA_test_PHN_Vdb_unc.mat','err_AOA_test_PHN_Vdb_unc');


    save('err_loc_test_PHN_Vdb_NoPHN.mat','err_loc_test_PHN_Vdb_NoPHN');
    save('err_ch_amp_test_PHN_Vdb_NoPHN.mat','err_ch_amp_test_PHN_Vdb_NoPHN');
    save('err_ch_test_PHN_Vdb_NoPHN.mat','err_ch_test_PHN_Vdb_NoPHN');
    save('err_AOA_test_PHN_Vdb_NoPHN.mat','err_AOA_test_PHN_Vdb_NoPHN');
    
    
    
    save('Wrong_num_PHN_Vdb.mat','Wrong_num_PHN_Vdb');
    save('small_h_num_PHN_Vdb.mat','small_h_num_PHN_Vdb');
    save('small_h_num_PHN_Vdb_NoPHN.mat','small_h_num_PHN_Vdb_NoPHN');
    
end


 



% % Debug
% MU_Pos;
% cord_est;
% norm(MU_Pos-cord_est);

% figure();stem(wrapToPi(PHN_est(:)));axis([0 P.N*P.Nr -0.5 0.5])
% figure();stem(wrapToPi(PHN(:)));axis([0 P.N*P.Nr -0.5 0.5])
% figure();stem(wrapToPi(PHN(:)-PHN_est(:)));axis([0 P.N*P.Nr -0.5 0.5])
%%
delete(gcp('nocreate'));
toc
























