%% OK, let's start with multipath coding!
clear;
close;
tic;
% parpool('local',8);               % Use parallel computting tool
%% Set up the system
P = setUpSystem();

% noise/power parameters
PHN_deg = 6;            % PHN std (in , typical value 0 deg ~ 10 deg)
var_x_dBm = 30;    % transmit power in dBm (  x = 10*log10(P/1mW), P = 1mW * 10^(x/10)  );

%% Position of the MU and Scatterers

% position of the user
row_index = randi([1 P.grid_side]);                         % pick one grid from x axis
col_index = randi([1 P.grid_side]);                         % pick one grid from y axis
MU_Pos_fixed = P.qs_mesh{row_index,col_index};              % pick out the rand grid
qs_mat = [P.qs{:}];                                         % 2 x 400 matrix containing the grid coordinate (left to right, up to down)
pos_index_tr = (row_index-1)*P.grid_side + col_index       % real index of the grid where the user is in
x_off_tr = unifrnd(-P.grid_length/2,P.grid_length/2)       % x axis off-grid
y_off_tr = unifrnd(-P.grid_length/2,P.grid_length/2)       % y axis off-grid
% x_off_tr = 0;      % x axis off-grid
% y_off_tr = 0;      % y axis off-grid
MU_pos_tr = MU_Pos_fixed + [x_off_tr;y_off_tr];             % true coordinate of the user

% position of the scatterers
scat_pos_tr = [[10;175] [-20;-225]];
% scat_pos_tr = [[-200;200] [200;200] [-200;-200] [200;-200]];

effect = [1 1 0 0;0 0 1 1].';
% effect = [1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1].';

K = size(scat_pos_tr, 2);

%% compute pathloss/receive sinr
pl_exp_LoS = 2;     % LoS pathloss exponent
pl_exp_NLoS = 3;    % NLoS pathloss exponent
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
y_err = y-Omega_tr*w_true -Phi_tr*v_true;

%% Start the algorithm of localization with PHN with NLoS and dynamic grid update
[cord_est,h_est,h_est_source] = loc_SRED_SPI_dyn(y,SNRs_LoS,SNRs_NLoS,PHN_deg,PHN_CovMat,P);

h_err = norm(h_est - h_tr)/norm(h_tr);
h_err_source = norm(h_est_source - h_tr)/norm(h_tr);
loc_err = norm(cord_est - MU_pos_tr);
fprintf("localization error is %f\n", loc_err)






























