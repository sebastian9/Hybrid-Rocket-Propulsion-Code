%% Combustion Chamber %%
function [result] = CombustionChamber(n_oxv_dot, n_oxl_dot, R_p, T_C, P_C, rho_C, M_ox_C, M_f_C, M_a_C, mi)

    m_dot_ox = -(n_oxl_dot - n_oxv_dot)*mi.MW_ox;
    m_dot_f = mi.R_p_dot(m_dot_ox)*(mi.L_g*pi*2*mi.R_p)*mi.rho_f;
    
    V_C = 2*pi*R_p^2*mi.L_g;
    V_C_dot = 2*pi*R_p*mi.L_g*mi.R_p_dot(m_dot_ox);

    M_C = M_ox_C + M_f_C + M_a_C;
%     OF = M_ox_C/M_f_C;
     
%     [~,~,k_C] = mi.NASACEA(OF, P_C);
    k_C = 1.33;
    lambda = sqrt(k_C*(2/(k_C+1))^((k_C+1)/(k_C-1)));
    m_dot_nz = lambda*mi.A_t*sqrt(P_C*rho_C);
    rho_C_dot = (m_dot_ox/V_C) + (mi.rho_f-rho_C)*(V_C_dot/V_C) ...
    	- (m_dot_nz/V_C);
    P_C_dot = (k_C-1)*mi.Q*(mi.rho_f)*(V_C_dot/V_C) ...
    	- P_C*k_C*(m_dot_nz/(V_C*sqrt(rho_C))) - (k_C-1)*mi.q/(V_C);

    result(1,1) = mi.R_p_dot(m_dot_ox);
    result(2,1) = P_C_dot;
    result(3,1) = T_C*((P_C_dot/P_C)-(rho_C_dot/rho_C));
    result(4,1) = rho_C_dot;
    result(5,1) = m_dot_ox - m_dot_nz*(M_ox_C/M_C);
    result(6,1) = m_dot_f - m_dot_nz*(M_f_C/M_C);
    result(7,1) = -m_dot_nz*(M_a_C/M_C);
end