function SigmaPointFilter(lidar,car)
    

    %Initialize vectors and matrices
    dt = 0.1;
    P0 = diag([2 5 1 pi/4 4 2])^2;
    Q = (diag([0.25 0.25 3 4*pi/180 0.1 0.1])^2)/dt;
    R = diag([2*pi/180 2*pi/180 0.1])^2;
    alpha = 10e-3;
    beta = 2;
    k = 1;
    x0 = [90 4.25 13 pi 5 2]';
    lidar_data = load('problem3data (1).mat');
    
    %calculate sigma points for propagation step
    function sig_vec = sigCreateProp(x_hat,P_hat)
        nx = size(x_hat,1);
        nv = 6;
        x_hat_aug = [x_hat; zeros(6,1)];
        P_aug = [P_hat zeros(6,6); zeros(6,6) Q];
        S = chol(P_aug)';
        lambda = alpha^2*(nx+nv+k)-(nx+nv);
        
        for i = 1:2*(nx+nv)+1
            if i == 1
                sig = x_hat_aug;
                sig_vec(:,i) = sig;
            elseif i <= nx+nv+1
                sig = x_hat_aug + sqrt(nx+nv+lambda)*S(:,i-1);
                sig_vec(:,i) = sig;
            else
                sig = x_hat_aug - sqrt(nx+nv+lambda)*S(:,i-1-(nx+nv));
                sig_vec(:,i) = sig;
            end
        end
    end

    %propagate sigma points 
    function eta_vec = sigProp(sig_vec)
        eta_vec = zeros(size(sig_vec,1),size(sig_vec,2));
        for i = 1:size(sig_vec,2)
            [eta] = ode45(@dyn_car,[0,dt],sig_vec(:,i));
            eta_vec(:,i) = eta.y(:,11);
        end  
    end

    %calculate x_bar and P_bar
    function [x_bar_aug, P_bar_aug] = CalcEtaSum(eta_vec)
        nx = 6;
        nv = 6;
        lambda = alpha^2*(nx+nv+k)-(nx+nv);
        
        
        W0_m = lambda/(nx+nv+lambda);
        W0_c = (lambda/(nx+nv+lambda))+1-alpha^2+beta;
        Wi_m = 1/(2*(nx+nv+lambda));
        Wi_c = 1/(2*(nx+nv+lambda));
        
        for iii = 1:2*(nx+nv)+1
            if iii == 1
                x_bar_temp1 = W0_m*eta_vec(:,iii);
                x_bar_vec(:,iii) = x_bar_temp1;
            else
                x_bar_temp2 = Wi_m*eta_vec(:,iii);
                x_bar_vec(:,iii) = x_bar_temp2;
            end
        end
        x_bar_aug = sum(x_bar_vec,2);

        for j = 1:2*(nx+nv)+1
            if j == 1
                P_temp = W0_c*(eta_vec(:,j)-x_bar_aug)*(eta_vec(:,j)-x_bar_aug)';
            else
                P_temp = P_temp + Wi_c*(eta_vec(:,j)-x_bar_aug)*(eta_vec(:,j)-x_bar_aug)';
            end
        end
        P_bar_aug = P_temp;
    end

    %extract lidar data ??????????????
    function [z_vec] = LidarToMeas(lidar)
        z_vec = zeros(3,size(lidar,1));
        for i = 1:size(lidar,1)
            min_r = min(lidar(i).z(:,1));
            max_b = max(lidar(i).z(:,2));
            min_b = min(lidar(i).z(:,2));
            
            z_temp = [min_b max_b min_r]';
            z_vec(:,i) = z_temp;
        end
    end

    %update with measurements
    function [gamma_vec,z_bar_i_vec] = sigCreateUpdate(x_bar_aug,P_bar_aug)
        x_bar = x_bar_aug(1:6,:);
        P_bar = P_bar_aug(1:6,1:6);
        
        x_bar_a = [x_bar; zeros(3,1)];
        P_bar_a = [P_bar zeros(6,3); zeros(3,6) R];
        
        nx = 6;
        nz = 3;
        Su = chol(P_bar_a)';
        
        lambda_u = alpha^2*(nx+k+3)-(nx+3);
        
        for i = 1:2*(nx+nz)+1
            if i == 1
                gamma = x_bar_a;
                gamma_vec(:,i) = gamma;
            elseif i <= (nx+nz)+1
                gamma = x_bar_a + sqrt(nx+nz+lambda_u)*Su(:,i-1);
                gamma_vec(:,i) = gamma;
            else
                gamma = x_bar_a - sqrt(nx+nz+lambda_u)*Su(:,i-1-(nx+nz));
                gamma_vec(:,i) = gamma;
            end
        end
        
        for ii = 1:2*(nx+nz)+1
            z_i = h_car(gamma_vec(:,ii));
            z_bar_i_vec(:,ii)=z_i;
        end
    end
    

    %calculate z_bar, Pzz, and Pxz
    function [z_bar,Pzz,Pxz] = Update(z_bar_i_vec,gamma_vec,x_bar_aug)
        nx = 6;
        nz = 3;
        lambda_u = alpha^2*(nx+k+nz)-(nx+nz);
        W0_m = lambda_u/(nx+nz+lambda_u);
        W0_c = lambda_u/(nx+nz+lambda_u) + 1 - alpha^2 + beta;
        Wi_m = 1/(2*(nx+nz+lambda_u));
        Wi_c = 1/(2*(nx+nz+lambda_u));
        x_bar = x_bar_aug(1:6);
        gamma_vals = gamma_vec(1:6,:);
        
        for i = 1:2*(nx+nz)+1
            if i == 1
                z_bar_temp1 = W0_m*z_bar_i_vec(:,i);
                z_bar_vec(:,i) = z_bar_temp1;
            else
                z_bar_temp2 = Wi_m*z_bar_i_vec(:,i);
                z_bar_vec(:,i) = z_bar_temp2;
            end
        end
        z_bar = sum(z_bar_vec,2);

        for j = 1:2*(nx+nz)+1
            if j == 1
                Pzz_temp = W0_c*(z_bar_i_vec(:,j)-z_bar)*(z_bar_i_vec(:,j)-z_bar)';
            else
                Pzz_temp = Pzz_temp + Wi_c*(z_bar_i_vec(:,j)-z_bar)*(z_bar_i_vec(:,j)-z_bar)';
            end
        end
        Pzz = Pzz_temp;
        
        for ii = 1:2*(nx+nz)+1
            if ii == 1
                Pxz_temp = W0_c*(gamma_vals(:,ii)-x_bar)*(z_bar_i_vec(:,ii)-z_bar)';
            else
                Pxz_temp = Pxz_temp + Wi_c*(gamma_vals(:,ii)-x_bar)*(z_bar_i_vec(:,ii)-z_bar)';
            end
        end
        Pxz = Pxz_temp;
    end

    %MMSE measurement update
    function [x_hat,Pxx_hat] = MMSE(z_bar,Pzz,Pxz,x_bar_aug,z_vec,ind,P_bar_aug)
        P_bar = P_bar_aug(1:6,1:6);
        x_bar = x_bar_aug(1:6);
        x_hat = x_bar + Pxz*inv(Pzz)*(z_vec(:,ind)-z_bar);
        Pxx_hat = P_bar - Pxz*inv(Pzz)*Pxz';
    end
    
    %%MAIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_vec = LidarToMeas(lidar);
    x_hat_vec = zeros(6,220);
    P_hat_vec = zeros(6,6,220);
    for jjj = 1:size(lidar,1)
        if jjj == 1
            sig_vec = sigCreateProp(x0,P0);
            eta_vec = sigProp(sig_vec);
            [x_bar_aug, P_bar_aug] = CalcEtaSum(eta_vec);
            [gamma_vec,z_bar_i_vec] = sigCreateUpdate(x_bar_aug,P_bar_aug);
            [z_bar,Pzz,Pxz] = Update(z_bar_i_vec,gamma_vec,x_bar_aug);
            [x_hat,Pxx_hat] = MMSE(z_bar,Pzz,Pxz,x_bar_aug,z_vec,jjj,P_bar_aug);
            x_hat_vec(:,jjj) = x_hat;
            P_hat_vec(:,:,jjj) = Pxx_hat;
        else
            sig_vec = sigCreateProp(x_hat,Pxx_hat);
            eta_vec = sigProp(sig_vec);
            [x_bar_aug, P_bar_aug] = CalcEtaSum(eta_vec);
            [gamma_vec,z_bar_i_vec] = sigCreateUpdate(x_bar_aug,P_bar_aug);
            [z_bar,Pzz,Pxz] = Update(z_bar_i_vec,gamma_vec,x_bar_aug);
            [x_hat,Pxx_hat] = MMSE(z_bar,Pzz,Pxz,x_bar_aug,z_vec,jjj,P_bar_aug);
            x_hat_vec(:,jjj) = x_hat;
            P_hat_vec(:,:,jjj) = Pxx_hat;
        end
    end
    
    for qq = 1:220
        plotcar(x_hat_vec(:,qq));
    end
    
    car_timehist = zeros(6,220);
    car_t_vec = zeros(1,220);
    for pp = 1:220
        car_x_val = car(pp).x;
        car_t_val = car(pp).t;
        car_timehist(:,pp) = car_x_val;
        car_t_vec(:,pp) = car_t_val;
    end
    
    final_estimate = x_hat_vec(:,220)
    actual_final = car_timehist(:,220)
    
    
end

    


        