function [F] = CV3 (AeAt, A_e, k_c, P_c, P_a, T_c, R, C_star, n_eff, theta_noz)
    
    A_t = A_e/AeAt;
    lamda = 0.5*(1+cos(theta_noz*pi/180));

    syms M_ei
    eqn = AeAt == (1/M_ei)*((1+(k_c-1)/2*M_ei^2)/(1+(k_c-1)/2))^((k_c+1)/(2*(k_c-1)));
    M_e_1st = double(vpasolve(eqn, M_ei, [0 1]));
    M_e_3rd = double(vpasolve(eqn, M_ei, [1 10]));
    P_e_1st = P_c*(1+(k_c-1)/2*M_e_1st^2)^(-k_c/(k_c-1));
    P_e_3rd = P_c*(1+(k_c-1)/2*M_e_3rd^2)^(-k_c/(k_c-1));
    P_1st = P_e_1st/P_c;
    P_3rd = P_e_3rd/P_c;
    M_1 = M_e_3rd;
    M_2 = ((M_1^2+2/(k_c-1))/(2*k_c/(k_c-1)*(M_1^2-1)))^0.5;
    P2P1 = (1+k_c*M_1^2)/(1+k_c*M_2^2);
    P_2nd = P2P1*P_3rd;
    
    
    if (P_a/P_c) >= P_1st
       
        M_e = (((P_c/P_a)^((k_c-1)/k_c)-1)/((k_c-1)/2))^0.5;
        T_e = T_c/(1+(k_c-1)/2*M_e^2);
        v_e = M_e*(k_c*R*T_e)^0.5;
        rho_e = P_a/R/T_e;
        m_dot_noz = rho_e*A_e*v_e/n_eff;
        
        F = lamda*n_eff*m_dot_noz*v_e;
        
    end
    
    if (P_a/P_c) < P_2nd
        
        M_e = M_e_3rd;
        P_e = P_c/(1+((k_c-1)/2)*M_e^2)^(k_c/(k_c-1));
        T_e = T_c/(1+(k_c-1)/2*M_e^2);
        v_e = M_e*(k_c*R*T_e)^0.5;
        m_dot_noz = P_c*A_t/n_eff/C_star;
        F = lamda*n_eff*m_dot_noz*v_e+(P_e-P_a)*A_e;
        
    end
    
    if (P_a/P_c) < P_1st && (P_a/P_c) >= P_2nd
        
        eqn = AeAt*P_e/P_c ==(1/M_ei)*((1+(k_c-1)/2*M_ei^2)/(1+(k_c-1)/2))^((k_c+1)/(2*(k_c-1)))*((1+((k_c-1)/2)*M_e^2)^(k_c/(k_c-1)))^-1;
        M_e = double(vpasolve(eqn, M_ei, [0 10]));
        P_o_e = P_a*(1+((k_c-1)/2)*M_e^2)^(k_c/(k_c-1));
        P_o_21 = P_o_e/P_c;
        syms M_1i
        eqn = P_o_21==(((k_c+1)/2*M_1i^2)/(1+(k_c-1)/2*M_1i^2))^(k_c/(k_c-1))*((2*k_c*M_1i^2/(k_c+1))-(k_c-1)/(k_c+1))^(1/(1-k_c));
        M_1 = double(vpasolve(eqn, M_1i, [0 10]));
        M_2 = ((M_1^2+2/(k_c-1))/(2*k_c/(k_c-1)*(M_1^2-1)))^0.5;
        T_1 = T_c/(1+(k_c-1)/2*M_1^2);
        T_2 = T_1*((1+(k_c-1)/2*M_1^2)/(1+(k_c-1)/2*M_2^2));
        T_o_2 = T_2*(1+(k_c-1)/2*M_2^2);
        T_e = T_o_2/(1+(k_c-1)/2*M_e^2);
        v_e = M_e*(k_c*R*T_e)^0.5;
        m_dot_noz = P_c*A_t/n_eff/C_star;
        F = lamda*n_eff*m_dot_noz*v_e;
        
    end
    
%     function [F] = CV3 ()
%         p_e = 100000;
%         p = 750000;
%         k = 1.4;
%         p_H = 100000;
%         sigma_c = 0.90;
%         A_e = 0.0158e-3;
%         p_er = p_e/p;
%         AtAe = (2/(k-1)*(2/(k+1))^((k+1)/(1-k))*p_er^(2/k)*(1-p_er^((k-1)/k)))^0.5;
%         AtAe = 0.397e-4/A_e;
%         a = (1-k)^(-1);
%         b = (1-k)/(1+k);
%         c = AtAe*((k+1)/2)^a;
% 
%         change = 1;
%         lamda = 2;
%         while round(change,30) > 0
%             lamda_plus = lamda-(lamda-c*(1+b*lamda^2)^a)/(1-2*a*b*c*lamda*(1+b*lamda^2)^(a-1));
%             change = ((lamda_plus - lamda)^2)^0.5;
%             lamda = lamda_plus
%         end
% 
%         F = A_e*p_H*(sigma_c*p/p_H*AtAe*k*(2/(k+1))^(1/(k+1))*lamda-1);
%     
%     end
end