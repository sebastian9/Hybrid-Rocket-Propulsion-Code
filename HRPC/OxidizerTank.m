%% Oxidizer Tank %%
function [result] = OxidizerTank(n_oxv, n_oxl, T_T, P_C, mi)
    
    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez
    
    P_T = (n_oxv + mi.n_spv) * mi.env.R_u * T_T / (mi.V_T - n_oxl*mi.V_mol_oxl(T_T));

    if n_oxl <= 0
        % Factors in the ODE system
        f1 = -mi.C_d_l * mi.N_inj * mi.A_inj * ...
            sqrt( 2* ( P_T - mi.P_losses - P_C ) / ( mi.MW_ox * 0.024789598e-3 ) );
        f2 = -mi.env.R_u * T_T;        
        f3 = mi.m_T * mi.c_P_T + ... 
            n_oxv * mi.CV_oxv(T_T) + mi.n_spv * mi.CV_spv(T_T);
        M = [1 0; f2 f3]; V = [f1; 0];
        X = linsolve(M,V);
        result = [X(1), 0, X(2)];
    else
        % Factors in the ODE system
        f1 = -mi.C_d_l * mi.N_inj * mi.A_inj * ...
            sqrt( 2* ( P_T - mi.P_losses - P_C ) / ( mi.MW_ox * mi.V_mol_oxl(T_T) ) );
        f2 = -mi.V_mol_oxl(T_T) * mi.P_crit_ox(T_T);
        f3 = ( mi.V_T - n_oxl * mi.V_mol_oxl(T_T) ) * mi.dP_crit_ox_dT(T_T) ...
            - n_oxl * mi.P_crit_ox(T_T) * mi.dV_mol_oxl_dT(T_T);
        f4 = mi.env.R_u * n_oxv;
        f5 = mi.env.R_u * T_T;
        f6 = mi.m_T * mi.c_P_T + n_oxl * CP_oxl + ... 
            n_oxv * mi.CV_oxv(T_T) + mi.n_spv * mi.CV_spv(T_T);
        f7 = mi.env.R_u * T_T - mi.deltaH_oxv(T_T);
        f8 = P_T * mi.V_mol_oxl(T_T);

        M = [ 1 1 0 ; f5 -f2 f4-f3 ; f7 f8 -f6 ]; % ODE System
        V = [ f1 ; 0 ; 0 ];

        [result] = linsolve(M,V);  % Isolated derivatives
        % n_oxv_dot, n_oxl_dot, T_T_dot
    end
end