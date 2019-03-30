%% State function definition
function state_dot = HRPC_function(t,state,env) 
    
    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
    R_p = state(4);
    rho_C = state(5);
    P_C = state(6);
    T_C = state(7);
   
    mi = env.Motor; % mi = motor_instance

    %% Oxidizer Tank %%
    
    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez
    
    P_T = (n_oxv + mi.n_spv) * env.R_u * T_T / (mi.V_T - n_oxl*mi.V_mol_oxl(T_T));

    % Factors in the ODE system
    f1 = -mi.C_d_l * mi.N_inj * mi.A_inj * ...
        sqrt( ( P_T - mi.P_losses - P_C ) / ( 2 * mi.MW_ox * mi.V_mol_oxl(T_T) ) );
    f2 = -mi.V_mol_oxl(T_T) * mi.P_crit_ox(T_T);
    f3 = ( mi.V_T - n_oxl * mi.V_mol_oxl(T_T) ) * mi.dP_crit_ox_dT(T_T) ...
        - n_oxl * mi.P_crit_ox(T_T) * mi.dV_mol_oxl_dT(T_T);
    f4 = env.R_u * n_oxv;
    f5 = env.R_u * T_T;
    f6 = mi.m_T * mi.c_P_T + n_oxl * CP_oxl + ... 
        n_oxv * mi.CV_oxv(T_T) + mi.n_spv * mi.CV_spv(T_T);
    f7 = env.R_u * T_T - mi.deltaH_oxv(T_T);
    f8 = P_T * mi.V_mol_oxl(T_T);

    M = [ 1 1 0 ; f5 -f2 f4-f3 ; f7 f8 -f6 ]; % ODE System
    V = [ f1 ; 0 ; 0 ];
    
    state_dot(1:3,1) = linsolve(M,V);  % Isolated derivatives
    
    %% Control Volume 2

    n_oxv_dot = state_dot(1);
    n_oxl_dot = state_dot(2);

    m_dot_ox = -(n_oxl_dot - n_oxv_dot)*mi.MW_ox;
    m_dot_f = mi.R_p_dot(m_dot_ox)*(mi.L_g*pi*2*mi.R_p)*mi.rho_f;
    
    V_C = 2*pi*R_p^2*mi.L_g;
    V_C_dot = 2*pi*R_p*mi.L_g*mi.R_p_dot(m_dot_ox);

    rho_ox_C = P_C/(env.R_u/T_C)
    V_ox_C = V_C*((rho-mi.rho_f)/(rho_ox_C-rho_f));
    V_f_C = V_C - V_ox_C;
    OF = (rho_ox_C*V_ox_C)/(mi.rho_f*V_f_C);
     
    [Q, k_C] = mi.NASACEA(OF);
    lambda = sqrt(k*(2/(k_C+1))^((k_C+1)/(k_C-1)));
    m_dot_nz = lambda*mi.A_t*sqrt(P_C*rho_C);
    rho_C_dot = (m_dot_ox/V_C) + (mi.rho_f-rho_c)*(V_C_dot/V_C) ...
    	+ (m_dot_nz/V_C);
    P_C_dot = (k_C-1)*Q*(mi.rho_f/P_C)*(V_C_dot/V_C) ...
    	- k_C*(m_dot_nz/(V_C*sqrt(rho))) - (k_C-1)*mi.q(T_C)/(P_C*V_C);

    state_dot(4) = mi.R_p_dot(m_dot_ox);
    state_dot(5) = rho_C_dot;
    state_dot(6) = P_C_dot;
    state_dot(7) = T_C*((P_C_dot/P_C)-(rho_C_dot/rho_C));

end