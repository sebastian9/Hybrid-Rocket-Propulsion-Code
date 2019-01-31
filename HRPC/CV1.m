syms A_inj C_d c_P_T CP_oxl CV_oxv CV_spv deltaH_oxv m_T MW_ox n_oxv_dot n_oxv n_oxl_dot n_oxl n_spv N_inj P_T P_losses P_C P_crit_ox dP_crit_ox_dT R_u T_T_dot T_T V_T V_mol_oxl dV_mol_oxl_dT

state_0 = [0 ; 0 ; 0];
tspan = 0:0.0005:5;

[t,state] = ode45(@CV1_function,tspan,state_0);

function state_dot = CV1_function(t,state)
    
    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);

    f1 = -C_d*N_inj*A_inj*sqrt((P_T-P_losses-P_C)/(2*MW_ox*V_mol_oxl));
    f2 = -V_mol_oxl*P_crit_ox;
    f3 = ((V_T - n_oxl*V_mol_oxl)*dP_crit_ox_dT - n_oxl*P_crit_ox*dV_mol_oxl_dT);
    f4 = R_u*n_oxv;
    f5 = R_u*T_T;
    f6 = m_T*c_P_T + n_oxl*CP_oxl + n_oxv*CV_oxv + n_spv*CV_spv;
    f7 = R_u*T_T - deltaH_oxv;
    f8 = P_T*V_mol_oxl;

    M = [ 1 1 0 ; f5 -f2 f4-f3 ; f7 f8 -f6 ];
    V = [ f1 ; 0 ; 0 ];

    state_dot = linsolve(M,V);

end

