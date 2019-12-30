function TruthModelSim()    
    clear all

    F     = [  0.81671934103521,  0.08791146849849;...
              -3.47061412053765,  0.70624978972000];     % for all k
    Gamma = [  0.00464254201630;...
               0.08791146849849];                        % for all k
    H     = [  2.00000000000000,  0.30000000000000];     % for all k
%
    Q     =    4.00000000000000;                         % for all k
    R     =    0.01000000000000;                         % for all k
%
    xbar0   = [  0.20000000000000;...
               -2.50000000000000];
    P0      = [  0.25000000000000,  0.08000000000000;...
                0.08000000000000,  0.50000000000000];
    xhat0   = [  0.20000000000000;...
               -2.50000000000000];


    function [zhist,xhist] = mcltisim(F,Gamma,H,Q,R,xbar0,P0,kmax)
        x_k = Init_state(P0,xbar0);
        xhist(1,1:2) = x_k;
  
        for i = 1:kmax
            z_k = H*x_k + Meas_Noise(R);
            zhist(i,1) = z_k;
            x_k = F*x_k + Gamma*Proc_Noise(Q);
            xhist(i+1,1:2) = x_k;
        end
    
        %local functions
        function x0 = Init_state(P0,xbar0)
            x0_ = randn(2,1);
            R_a = chol(P0);
            x0 = transpose(R_a)*x0_+xbar0;
        end
    
        function v_k = Proc_Noise(Q)
            v_k_ = randn(1);
            R_a = chol(Q);
            v_k = transpose(R_a)*v_k_;
        end
    
        function w_k = Meas_Noise(R)
            w_k_ = randn(1);
            R_a = chol(R);
            w_k = transpose(R_a)*w_k_;
        end
    end

    function [x_hat_vec,p_hat_vec] = Kalman(F,Gamma,H,Q,R,P0,zhist,xhat0)
        x_bar = F*xhat0;
        p_bar = F*P0*F' + Gamma*Q*Gamma';
        x_hat_vec(1,1:2) = xhat0;
        
        for i = 1:size(zhist)
            U = zhist(i) - H*x_bar;
            S = H*p_bar*H' + R;
            W = p_bar*H'*inv(S);
            x_hat = x_bar + W*U;
            p_hat = p_bar - W*S*W';
            
            x_hat_vec(i,1:2) = x_hat;
            p_hat_vec(1:2,1:2,i) = p_hat;
            
            x_bar = F*x_hat;
            p_bar = F*p_hat*F' + Gamma*Q*Gamma';
            
        end
        x_hat_vec(36,1:2) = x_bar;
    end

    function [x_err_tens,p_hat_vec] = DoMonteCarlo(N,F,Gamma,H,Q,R,xbar0,P0,xhat0)
        for k = 1:N
            [zhist,xhist] = mcltisim(F,Gamma,H,Q,R,xbar0,P0,35);
            [x_hat_vec,p_hat_vec] = Kalman(F,Gamma,H,Q,R,P0,zhist,xhat0);
            x_err = xhist - x_hat_vec;
            x_err_tens(1:36,1:2,k) = x_err;
        end
    end

    function P_x_err_all = calcP_x(x_err_tens,N,row)
        for k = 1:N
            P_x_err = x_err_tens(row,1:2,k)'*x_err_tens(row,1:2,k);
            P_x_err_all(1:2,1:2,k) = P_x_err;
        end
        P_x_err = P_x_err;
    end

    [x_err_tens,p_hat_vec] = DoMonteCarlo(50,F,Gamma,H,Q,R,xbar0,P0,xhat0);
    E_x_err_all = sum(x_err_tens,3)/50;
    P_x_err_all_10 = calcP_x(x_err_tens,50,10);
    P_x_err_all_35 = calcP_x(x_err_tens,50,35);
    
    E_x_err_10 = E_x_err_all(10,1:2)
    E_x_err_35 = E_x_err_all(35,1:2)
    
    P_x_err_10 = sum(P_x_err_all_10,3)/50
    P_x_err_35 = sum(P_x_err_all_35,3)/50
    
    P_x_10 = p_hat_vec(:,:,10)
    P_x_35 = p_hat_vec(:,:,35)
    
    

    

end
    
        