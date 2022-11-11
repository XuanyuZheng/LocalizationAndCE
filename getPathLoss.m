function rhos = getPathLoss(MU_Pos,P,expn)
%getPathLoss get Path Loss Coefficients of N BSs
%   MU_Pos is the user coordinate, P is the system parameters
    PLs = zeros(P.N,1);
    PL0 = 20*log10(4*pi*P.fc/P.c);
    for n = 1:P.N
        distance_MU_BS = norm(P.BS_Pos(:,n)-MU_Pos);
        PLs(n) = PL0 + 10*expn*log10(distance_MU_BS);
    end
    rhos = 1./sqrt(10.^(PLs./10));                      % PL coefficients in N BSs
end

