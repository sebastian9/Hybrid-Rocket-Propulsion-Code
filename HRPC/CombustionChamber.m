%% Combustion Chamber %%
function [result] = CombustionChamber(n_oxv_dot, n_oxl_dot, R_p, T_C, P_C, rho_C, M_ox_C, M_f_C, M_a_C, mi)

    if n_oxl_dot ~= 0 % Liquid oxidizer is not depleted
        % The minus sign is due to the derivative being referred to the
        % oxidizer tank rather than to the combustion chamber; the oxidizer
        % tank losses matter whilst the combution chamber gains it. The
        % substraction comes from the fact that the gaseous oxidizer comes
        % from liquid oxidizer, but it is not fed to the combustion chamber
        % until the liquid is depleted.
        m_dot_ox = -(n_oxl_dot - n_oxv_dot)*mi.MW_ox; % Oxidizer mass flow rate [kg/s]
    else
        m_dot_ox = -n_oxv_dot*mi.MW_ox; % Oxidizer mass flow rate [kg/s]
    end
     
    V_C = pi*R_p^2*mi.L_g; % Combustion chamber volume [m^3]
    % Combustion chamber volume rate of change [m^3/s] A*dR/dt
    V_C_dot = 2*pi*R_p*mi.L_g*mi.R_p_dot(m_dot_ox,R_p,P_C); 
    m_dot_f = V_C_dot*mi.rho_f; % Fuel mass flow rate [kg/s] Rho*A*dR/dt

    M_C = abs(M_ox_C) + abs(M_f_C) + abs(M_a_C); % Combustion chamber mass [kg]

    k_C = 1.4; % Ratio of specific heats
    % Nozzle Control Volume, [kg/s], [N]
    [m_dot_nz,result(8,1)] = Nozzle(mi, k_C, P_C, rho_C);
    % Combustion chamber density rate of change [kg.m^-3.s^-1]
    % Chelaru, 2011, eq. 15
    rho_C_dot = (m_dot_ox/V_C) + (mi.rho_f-rho_C)*(V_C_dot/V_C) ...
    	- (m_dot_nz/V_C);
    % Combustion chamber pressure rate of change [Pa.s^-1] 
    % Derived from Chelaru, 2011, eq. 26 & 27
    P_C_dot = (k_C-1)*mi.Q*(mi.rho_f)*(V_C_dot/V_C) ...
    	- P_C*k_C*(m_dot_nz/(V_C*rho_C)) - (k_C-1)*mi.q/(V_C);

    result(1,1) = mi.R_p_dot(m_dot_ox,R_p,P_C);
    result(2,1) = P_C_dot;
    result(3,1) = T_C*((P_C_dot/P_C)-(rho_C_dot/rho_C)); % Chelaru, 2011, eq. 27 ! Replace for implementation of eq. 29
    result(4,1) = rho_C_dot;
    result(5,1) = 0;%m_dot_ox - m_dot_nz*(M_ox_C/M_C);
    result(6,1) = 0;%m_dot_f - m_dot_nz*(M_f_C/M_C);
    result(7,1) = 0;%-m_dot_nz*(M_a_C/M_C);
end