function pl_tilde = getAntennaCoordinate(type,M,L,lambda)

    % antenna coordinate relative to antenna center
    % same geometry for all BSs
    pl_tilde=zeros(M,2,L); % same for all BSs

    % BS_Ori=rand(L,1)*2*pi; %orientation of BS antennas
    BS_Ori = 0;
    %BS_Ori=pi/2*ones(L,1);
    for l=1:L
        switch type
            case 'ULA'
                arrayULA = phased.ULA('NumElements',M,'ElementSpacing',0.5*lambda);
                elementPos = getElementPosition(arrayULA); % along y axis
                pl_tilde(:,:,l) = elementPos(1:2,:)';
            case 'UCA'  % distance of two adjacent antennas set to r
                r = lambda/(2*sin(pi/M))/2;   % half lambda
        %        r=lambda/(2*sin(pi/M))/2;
        %        r=lambda;
                arrayUCA = phased.UCA('NumElements',M,'Radius',r);
                elementPos = getElementPosition(arrayUCA);
                pl_tilde(:,:,l) = elementPos(1:2,:)';
            case 'URA'   
                arrayURA = phased.URA('Size',[16,M/16],'ElementSpacing',0.5*lambda,'ArrayNormal','z'); 
                elementPos = getElementPosition(arrayURA);
                pl_tilde = elementPos(1:2,:)';
            case 'Random'%ULA
        %    pl_tilde=zeros(M,2); % antenna coordinate relative to antenna center
        for m=1:M/2
            pl_tilde(m,:,l) = [(lambda/4+lambda/2*(M/2-m))*cos(BS_Ori(l)),(lambda/4+lambda/2*(M/2-m))*sin(BS_Ori(l))];
        end
        
        for m=M/2+1:M
           pl_tilde(m,:,l) = [(lambda/4+lambda/2*(m-M/2-1))*cos(BS_Ori(l)+pi),(lambda/4+lambda/2*(m-M/2-1))*sin(BS_Ori(l)+pi)];  
        end 


        end
    end

end