function [a,theta] = steering_vec(dest_coord,source_coord,ant_coord,Nr,lambda_c)
    a = zeros(Nr,1);                                    %steering vector
    antenna_pos = ant_coord;
    
    theta = atan((source_coord(2)-dest_coord(2))/(source_coord(1)-dest_coord(1)));              % AOA
    if source_coord(1) < dest_coord(1)
       theta = theta + pi; 
    end
    for i = 1:Nr
        this_ant_pos = antenna_pos(:,i);
        a(i) = exp((2*pi*1i/lambda_c)*(this_ant_pos'*[cos(theta);sin(theta)]));
    end
end

