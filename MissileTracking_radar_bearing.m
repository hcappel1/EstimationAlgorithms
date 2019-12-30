function MissileTracking_radar_bearing(rhoahist,rhobhist,thetaahist,thetabhist,thist)
    %%initial guess
    l_a = 4.1e5;
    l_b = 4.4e5;
    
    for k = 1:12
        cos_theta = -0.5*(rhoahist(k)^2-rhobhist(k)^2-(l_b-l_a)^2)/(rhobhist(k)*(l_b-l_a));
        theta = acos(cos_theta);
        y_1_k = l_b - rhobhist(k)*cos(theta);
        y_2_k = rhobhist(k)*sin(theta);
   
        y_1(k) = y_1_k;
        y_2(k) = y_2_k;
    end
    y_1 = transpose(y_1);
    y_2 = transpose(y_2);
    
    v_1_0 = (y_1(12) - y_1(2))/(thist(12)-thist(2));
    v_2_0 = (y_2(12) - y_2(2) + 0.5*9.8*(thist(12)^2 - thist(2)^2))/(thist(12)-thist(2));
    y_1_0 = y_1(2)-v_1_0*thist(2);
    y_2_0 = y_2(2)-v_2_0*thist(2);
   
    %x_g = [y_1_0;v_1_0;y_2_0;v_2_0];
    x_g = 1.0e+03 *[
    1.5062;
    0.8997;
    1.6849;
    1.4998;]

    %%
    z = data(rhoahist,rhobhist,thetaahist,thetabhist);
    R = calcR();
    H = calc_H(x_g,thist)
    rank(H)
    %size(H);
    %size(R);
    %[x_hat,J_cost,P_xx] = GaussNewton(x_g,z,thist)
    
    %define local functions
    function R = calcR()
        R_small = [100 0 0 0; 0 900 0 0; 0 0 (.01)^2 0; 0 0 0 (.03)^2];
        l = 0;
        for n = 1:28
            R(n+3*l:n+3*l+3,n+3*l:n+3*l+3) = R_small;
            l = l + 1;
        end
    end

    function z = data(rhoahist,rhobhist,thetaahist,thetabhist)
        n = 0;
        for i = 1:28
            z(i+3*n) = rhoahist(i);
            z(i+3*n+1) = rhobhist(i);
            z(i+3*n+2) = thetaahist(i);
            z(i+3*n+3) = thetabhist(i);
            n = n+1;
        end
        z = transpose(z);
    end

    function [x_hat,J_cost,P_xx] = GaussNewton(x_g,z,t)
        R = calcR();
        
        delta_x = [1;1;1;1];
        while norm(delta_x) > 1e-8
            J = 0;
            norm_d_x = norm(delta_x);
            n = 0;
            H = calc_H(x_g,t);
            for i = 1:28
                h_temp = calc_h(x_g,t(i));
                q = [z(i+3*n); z(i+3*n+1); z(i+3*n+2); z(i+3*n+3)] - h_temp;
                d(i+3*n) = q(1);
                d(i+3*n+1) = q(2);
                d(i+3*n+2) = q(3);
                d(i+3*n+3) = q(4);
                n = n+1;
            end
            l = transpose(d);
            delta_x = inv(transpose(H)*H)*transpose(H)*l;
            x_g = x_g + delta_x;
            
        end
        x_hat = x_g;
        J_cost = transpose(l)*inv(R)*l;
        P_xx = inv(transpose(H)*inv(R)*H);
    end
        
    function h = calc_h(x_g,t_i)
        l_a = 4.1e5;
        l_b = 4.4e5;
        
        d_y_1_a = l_a-x_g(1)-t_i*x_g(2);
        d_y_1_b = l_b-x_g(1)-t_i*x_g(2);
        d_y_2 = x_g(3)+t_i*x_g(4)-0.5*9.81*t_i^2;
        
        h(1) = sqrt((d_y_1_a)^2+(d_y_2)^2);
        h(2) = sqrt((d_y_1_b)^2+(d_y_2)^2);
        h(3) = atan2(d_y_2, d_y_1_a);
        h(4) = atan2(d_y_2, d_y_1_b);
        h = transpose(h);
    end

    function H = calc_H(x_g,t)
        l_a = 4.1e5;
        l_b = 4.4e5;
        
        H = zeros(112,4);
        n = 0;
        for i = 1:28
            d_y_1_a = l_a-x_g(1)-t(i)*x_g(2);
            d_y_1_b = l_b-x_g(1)-t(i)*x_g(2);
            d_y_2 = x_g(3)+t(i)*x_g(4)-0.5*9.81*t(i)^2;
            
            h_local = calc_h(x_g,t(i));
            
            H_11 = -d_y_1_a/h_local(1);
            H_12 = d_y_1_a/h_local(1)*(-t(i));
            H_13 = d_y_2/h_local(1);
            H_14 = d_y_2/h_local(1)*(t(i));
            
            H_21 = -d_y_1_b/h_local(2);
            H_22 = d_y_1_b/h_local(2)*(-t(i));
            H_23 = d_y_2/h_local(2);
            H_24 = d_y_2/h_local(2)*(t(i));
            
            H_31 = d_y_2/h_local(1)^2;
            H_32 = d_y_2/h_local(1)^2*t(i);
            H_33 = d_y_2/h_local(1)^2;
            H_34 = d_y_2/h_local(1)^2*t(i);
            
            H_41 = d_y_2/h_local(2)^2;
            H_42 = d_y_2/h_local(2)^2*t(i);
            H_43 = d_y_2/h_local(2)^2;
            H_44 = d_y_2/h_local(2)^2*t(i);
            
            H(i+3*n,1) = H_11;
            H(i+3*n,2) = H_12;
            H(i+3*n,3) = H_13;
            H(i+3*n,4) = H_14;
            
            H(i+3*n+1,1) = H_21;
            H(i+3*n+1,2) = H_22;
            H(i+3*n+1,3) = H_23;
            H(i+3*n+1,4) = H_24;
            
            H(i+3*n+2,1) = H_31;
            H(i+3*n+2,2) = H_32;
            H(i+3*n+2,3) = H_33;
            H(i+3*n+2,4) = H_34;
            
            H(i+3*n+3,1) = H_41;
            H(i+3*n+3,2) = H_42;
            H(i+3*n+3,3) = H_43;
            H(i+3*n+3,4) = H_44;
            
            n = n+1;
        end
    end
end