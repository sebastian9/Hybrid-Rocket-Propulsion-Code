%% State function definition
function state_dot = HRPC_function(t,state,env)  % mi = motor_instance
    
    % Curve fitted combustion chamber pressure [Pa]: 
    % P_C = -2924.42*t^6 + 46778.07*t^5 - 285170.63*t^4 + 813545.02*t^3 - ... 
    %    1050701.53*t^2 + 400465.85*t + 1175466.2;
    P_C = 500000; % Pa

    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
    
    mi = env.Motor;
    
    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez
    
    P_T = (n_oxv + mi.n_spv) * env.R_u * T_T / (mi.V_T-n_oxl * mi.V_mol_oxl(T_T));

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

    state_dot = linsolve(M,V); % Isolated derivatives 
end