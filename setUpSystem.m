function P = setUpSystem()
% System parameters can be changed inside this function only

    P = struct;
    
    % specifications on the BSs
    P.Nr = 30;                                          % # of anttenas of one BS
    P.fc = 30e9;                                        % carrier frequency in Hz
    P.c = 299792458;                                    % speed of light in m/s
    P.lambda_c = P.c/P.fc;                              % carrier wavelength in m
%     P.std_n = sqrt(P.var_n);   % noise standard deviation (no need if use receive SNR?)
    P.d = P.lambda_c/2;                                 % anteanna spacing in m
    P.x = 1;                                            % unit transmitted symbol, C^(1,1)
    P.antenna_pos = ...
    getAntennaCoordinate('UCA',P.Nr,1,P.lambda_c)';   % coordinates of the antennas relative to array center

    % noise power
    P.BW = 3.84e6;                                      % Band width in Hz
    P.Kb = 1.38e-23;                                    % Boltzmann constant
    P.T = 300;                                          % Temperature in Kelvin
    P.var_n = P.Kb*P.T*P.BW;                            % AWGN power in watts / variance of AWGN
    
    % specification of the map
    P.total_length = 100;                               % side length of the square map in m
    P.grid_side = 20;                                   % grid number on one side
    P.grid_length = P.total_length/P.grid_side;         % length of one grid
    P.dis_tm = 1.5;                                     % times of map total length for BS distance
    P.cornerVal = P.total_length*P.dis_tm;              % BS side distance
    P.BS_Pos = [-P.cornerVal P.cornerVal -P.cornerVal P.cornerVal;  % Four BS positions are fixed
           P.cornerVal P.cornerVal -P.cornerVal -P.cornerVal];  % BSs are at the four corners
    
    % coordinates of the fixed grids
    P.qs_mesh = cell(P.grid_side,P.grid_side);
    P.gap = P.grid_length/2;                                                % half grid length
    first_center = [-P.total_length/2+P.gap; P.total_length/2-P.gap];       % coordinate of the center of the first gird
    for i = 1:P.grid_side
       for j = 1:P.grid_side
           P.qs_mesh{i,j} = first_center + [(j-1)*2*P.gap;-(i-1)*2*P.gap];
       end
    end
    P.qs = vec((P.qs_mesh)');
    
    % coordinate of the fixed
    P.M = P.Nr;             % number of phi grid
    phi_intvl = 2*pi/P.M;   % grid size of phi
    phi_head = -pi/2;       % head in rad
    phi_tail = 3*pi/2-phi_intvl;      % tail in rad
    P.phi = (phi_head:phi_intvl:phi_tail).';    % grid edges
    P.phi_mesh = P.phi + phi_intvl/2;   % coordinate of the center of the girds
    
    % Others
    P.PL_exp = 2;                                       % Path Loss exponent
    P.Q = P.grid_side^2;                                % # of grids
    [~,P.N] = size(P.BS_Pos);                           % # of BSs

end

