function [m_dot_nz, F] = CV3 (mi, k_C, P_C, rho_C)
    
    A_t = mi.A_t;
    A_e = mi.A_e;
    n_eff = mi.n_eff;
    P_a = mi.env.P;
    theta_noz = mi.theta_noz;
    AtAe = A_t/A_e;
    sigma_c = 0.5*(1+cos(theta_noz*pi/180));

    P_1st = (2/(k_C+1))^(k_C/k_C-1); % Choked Flow Criteria
    
    if (P_a/P_C) >= P_1st % Nozzle is not choked
        v_e = sqrt(2*rho_C*(P_C-P_a));
        m_dot_nz = sigma_c*rho_C*A_e*v_e/n_eff;
        F = m_dot_nz*v_e;      
    else % Nozzle is choked
        m_dot_nz = sqrt(k_C*(2/(k_C+1))^((k_C+1)/(k_C-1)))*mi.A_t*sqrt(P_C*rho_C);
        a = (1-k_C)^(-1);
        b = (1-k_C)/(1+k_C);
        c = AtAe*((k_C+1)/2)^a;
        change = 1;
        lambda = 2;
        while change > 0.001
            lambda_plus = lambda-(lambda-c*(1+b*lambda^2)^a)/(1-2*a*b*c*lambda*(1+b*lambda^2)^(a-1));
            change = abs(lambda_plus - lambda);
            lambda = lambda_plus;
        end
        F = A_e*P_a*(sigma_c*P_C/P_a*AtAe*k_C*(2/(k_C+1))^(1/(k_C+1))*lambda-1);
    end
end