function a = steering_vec_NLoS(theta,ant_coord,Nr,lambda_c)
    a = zeros(Nr,1);                                    %steering vector
    antenna_pos = ant_coord;
    
    for i = 1:Nr
        this_ant_pos = antenna_pos(:,i);
        a(i) = exp((2*pi*1i/lambda_c)*(this_ant_pos'*[cos(theta);sin(theta)]));
    end
end

