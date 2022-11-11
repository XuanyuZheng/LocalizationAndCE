function rhos = getPathLossNLoS(MU_Pos,scat_Pos, P, expn)
%getPathLoss get Path Loss Coefficients of N BSs
%   MU_Pos is the user coordinate, P is the system parameters
    PLs = zeros(P.N,1);
    PL0 = 20*log10(4*pi*P.fc/P.c);
    for n = 1:P.N
%         distance_MU_scat = norm(scat_Pos-MU_Pos);
%         distance_scat_BS = norm(P.BS_Pos(:,n)-scat_Pos);
%         dist_overall = distance_MU_scat + distance_scat_BS;
        dist_overall = norm(P.BS_Pos(:,n)-MU_Pos);
        PLs(n) = PL0 + 10*expn*log10(dist_overall);
    end
    rhos = 1./sqrt(10.^(PLs./10));                      % PL coefficients in N BSs
end

