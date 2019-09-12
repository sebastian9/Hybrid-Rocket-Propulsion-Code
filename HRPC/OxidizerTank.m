%% Oxidizer Tank %%
function [result] = OxidizerTank(n_oxv, n_oxl, T_T, P_C, mi)

    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez 2009, Appendix H
    % Tank pressure, Ideal gas law of partial pressures
    P_T = (n_oxv + mi.n_spv) * mi.env.R_u * T_T / (mi.V_T - n_oxl*mi.V_mol_oxl(T_T));

    if n_oxl > 0
        % Factors in the ODE system. Fernanez, 2009 A.16
        f1 = -mi.C_d_l * mi.N_inj * mi.A_inj * ...
            sqrt( 2* ( P_T - mi.P_losses - P_C ) / ( mi.MW_ox * mi.V_mol_oxl(T_T) ) );
        f2 = -mi.V_mol_oxl(T_T) * mi.P_crit_ox(T_T);
        f3 = ( mi.V_T - n_oxl * mi.V_mol_oxl(T_T) ) * mi.dP_crit_ox_dT(T_T);
        f4 = mi.env.R_u * n_oxv;
        f5 = mi.env.R_u * T_T;
        f6 = mi.m_T * mi.c_P_T + n_oxl * CP_oxl + ...
            n_oxv * mi.CV_oxv(T_T) + mi.n_spv * mi.CV_spv(T_T);
        f7 = mi.env.R_u * T_T - mi.deltaH_oxv(T_T);
        f8 = P_T * mi.V_mol_oxl(T_T);

        % ODE System, Fern�ndez, 2009 A.16
        M = [ 1 1 0 ; f5 -f2 f4-f3 ; f7 f8 -f6 ];
        V = [ f1 ; 0 ; 0 ];
        % Isolated derivatives
        [result] = linsolve(M,V);
        % n_oxv_dot, n_oxl_dot, T_T_dot
    else % Liquid oxidizer depleted
        % Factors in the ODE system, G�nevieve 2013, 3.17 & 3.18
        f1 = -mi.C_d_l * mi.N_inj * mi.A_inj * ...
            sqrt( 2* ( P_T - mi.P_losses - P_C ) / ( mi.MW_ox * mi.env.R_u*T_T/P_T ) );
        f2 = -mi.env.R_u * T_T;
        f3 = mi.m_T * mi.c_P_T + ...
            n_oxv * mi.CV_oxv(T_T) + mi.n_spv * mi.CV_spv(T_T);
        M = [1 0; f2 f3]; V = [f1; 0];
        X = linsolve(M,V);
        % Isolated derivatives
        result = [X(1), 0, X(2)];
    end
end
