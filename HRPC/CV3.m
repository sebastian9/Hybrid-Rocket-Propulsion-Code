function [m_dot_nz, F] = CV3(mi, k_C, P_C, rho_C)
    AtAe = mi.A_t/mi.A_e;
    sigma_c = 0.9;%0.5*(1+cos(mi.theta_noz*pi/180));

    P_1st = (2/(k_C+1))^(k_C/(k_C-1)); % Choked Flow Criteria, Isentropic Flow Equations
    
    if mi.env.P/P_C >= P_1st % Nozzle is not choked
        if P_C-mi.env.P > 0
            % Bernoulli equation, exit static pressure is equal to
            % atmospheric
            v_e = sqrt(2*rho_C*(P_C-mi.env.P)); 
        else
            % Not backflow allowed in the model
            v_e = 0;
        end
        % Continuity equation
        m_dot_nz = sigma_c*rho_C*mi.A_e*v_e;
        % Thrust equation
        F = m_dot_nz*v_e;
    else % Nozzle is choked
        % Isentropic flow equations, Chelaru, 2011, eq. 12
        m_dot_nz = sqrt(k_C*(2/(k_C+1))^((k_C+1)/(k_C-1)))*mi.A_t*sqrt(P_C*rho_C);
        if m_dot_nz < 0; m_dot_nz = 0; end % No backflow allowed
        % Velocity ratio trascendental equation iterative solution, 
        % Chelaru, 2011, eq. 35
        lambda = fzero(@(lambda) VelocityFunction(lambda,k_C, AtAe),2); % Velocity ratio
        % Thrust, Chelaru, 2011, eq. 42
        F = mi.A_e*mi.env.P*(sigma_c*(P_C/mi.env.P)*AtAe*k_C*(2/(k_C+1))^(1/(k_C+1))*lambda-1);
        if F < 0; F = 0; end
    end
end