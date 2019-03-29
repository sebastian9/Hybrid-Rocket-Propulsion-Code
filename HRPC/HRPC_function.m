%% State function definition
function state_dot = HRPC_function(t,state,env) 
    
    mi = env.Motor; % mi = motor_instance
    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
%     % Curve fitted combustion chamber pressure [Pa]: 
%     P_C = -2924.42*t^6 + 46778.07*t^5 - 285170.63*t^4 + 813545.02*t^3 - ... 
%        1050701.53*t^2 + 400465.85*t + 1175466.2;
    R_p = state(4); % m 
    P_C = state(5); % Pa
    M_ox_C = state(6);
    M_f_C = state(7);
    %P_C_2_int = state(8);
    
    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez
    
    P_T = (n_oxv + mi.n_spv) * env.R_u * T_T / (mi.V_T - n_oxl*mi.V_mol_oxl(T_T));

    % Factors in the ODE system
    f1 = -mi.C_d_l * mi.N_inj * mi.A_inj * ...
        sqrt( 2 * ( P_T - mi.P_losses - P_C ) / ( mi.MW_ox * mi.V_mol_oxl(T_T) ) );
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
    
    n_oxv_dot = state_dot(1);
    n_oxl_dot = state_dot(2);
    
    m_dot_ox = -(n_oxl_dot - n_oxv_dot)*mi.MW_ox;
    m_dot_f = mi.R_p_dot(m_dot_ox,R_p)*(mi.L_g*pi*2*R_p)*mi.rho_f;
    
    state_dot(4) = mi.R_p_dot(m_dot_ox,R_p);
    
    V_C = mi.L_g*pi*R_p^2;
    dV_Cdt = 2*pi*R_p*mi.L_g*mi.R_p_dot(m_dot_ox,R_p);
    
    m_dot_n = 5 -.05*t;
    
    OF = (M_ox_C + m_dot_ox*0.001)/(M_f_C + m_dot_f*0.001);
    state_dot(6) = m_dot_ox - m_dot_n / ( 1 + OF^-1);
    state_dot(7) = m_dot_f - m_dot_n / ( 1 + OF);
    
    [T_C, CP_C, k_C] = mi.NASACEA(OF, P_C);
    
    state_dot(5) = ( 1000*CP_C*T_C*(m_dot_ox + m_dot_f - m_dot_n)*(k_C-1)/k_C - P_C*k_C*dV_Cdt ) / V_C;
    
    %state_dot(8) = (M_ox_C/mi.MW_ox + M_f_C/450.86)*env.R_u*T_C/V_C;
    
%     [state_dot(6),state_dot(7)] = fsolve(@mi.MassRate(M_dot,M_ox_C,M_f_C,m_dot_ox,m_dot_f,m_dot_n),[0,0]);
    
end