function [m_dot_nz, F] = Nozzle(mi, k_C, P_1, T_1)
    R = 259.8;
    C_F = (2*k_C/(k_C-1)*(2/(k_C+1))^((k_C+1)/(k_C-1))*(1-(mi.env.P/P_1)^((k_C-1)/k_C)))^0.5;
    F = C_F*mi.A_t*P_1;
    m_dot_nz = mi.A_t*P_1*k_C*((2/(k_C+1))^((k_C+1)/(k_C-1)))^0.5/(k_C*R*T_1)^0.5;
end
