%% State function definition %%
function state_dot = HRPC_function(t,state,env) 
    %% State Variables %% 
    
    n_oxv = state(1); 
    n_oxl = state(2);
    T_T = state(3);
    R_p = state(4);
    P_C = state(5);
    T_C = state(6);
    rho_C = state(7);
    M_ox_C = state(8);
    M_f_C = state(9);
    M_a_C = state(10);
   
    mi = env.Motor; % mi = motor_instance

    %% Oxidizer Tank %%
%     P_C = 50*10000; %Pa
    state_dot(1:3,1) = OxidizerTank(n_oxv, n_oxl, T_T, P_C, mi);
    %% Combustion Chamber %%

    n_oxv_dot = state_dot(1);
    n_oxl_dot = state_dot(2);

    m_dot_ox = -(n_oxl_dot - n_oxv_dot)*mi.MW_ox;
    m_dot_f = mi.R_p_dot(m_dot_ox)*(mi.L_g*pi*2*mi.R_p)*mi.rho_f;
    
    V_C = 2*pi*R_p^2*mi.L_g;
    V_C_dot = 2*pi*R_p*mi.L_g*mi.R_p_dot(m_dot_ox);

    M_C = M_ox_C + M_f_C + M_a_C;
    OF = M_ox_C/M_f_C;
     
%     [~,~,k_C] = mi.NASACEA(OF, P_C);
    k_C = 1.33;
    lambda = sqrt(k_C*(2/(k_C+1))^((k_C+1)/(k_C-1)));
    m_dot_nz = lambda*mi.A_t*sqrt(P_C*rho_C);
    rho_C_dot = (m_dot_ox/V_C) + (mi.rho_f-rho_C)*(V_C_dot/V_C) ...
    	+ (m_dot_nz/V_C);
    P_C_dot = (k_C-1)*mi.Q*(mi.rho_f)*(V_C_dot/V_C) ...
    	- P_C*k_C*(m_dot_nz/(V_C*sqrt(rho_C))) - (k_C-1)*mi.q/(V_C);

    state_dot(4) = mi.R_p_dot(m_dot_ox);
    state_dot(5) = P_C_dot;
    state_dot(6) = T_C*((P_C_dot/P_C)-(rho_C_dot/rho_C));
    state_dot(7) = rho_C_dot;
    state_dot(8) = m_dot_ox - m_dot_nz*(M_ox_C/M_C);
    state_dot(9) = m_dot_f - m_dot_nz*(M_f_C/M_C);
    state_dot(10) = -m_dot_nz*(M_a_C/M_C);

end