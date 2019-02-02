% Consideradas: A_inj C_d c_P_T CP_oxl CV_oxv CV_spv m_T MW_ox N_inj
% dP_crit_ox(T_T)_dT P_crit_ox(T_T) P_losses R_u T_T deltaH_oxv(T_T)
% No consideradas: n_oxv n_oxl n_spv P_T P_CV_T V_mol_oxl dV_mol_oxl_dT
clear

% %% Primal Variables
% 
% C_d_l = 0.8;
% % C_d_g = 0.5;
% C_d = C_d_l; 
% D_inj = 0.002; %m
% P_losses = 5; % bar
% N_inj = 17;
% m_T = 6.4882; % 'kg'
% m_oxt = 30; % 'kg'
% m_spv = 0.004; % 'kg'
% MW_ox = 44.013; % 'kg/kmol'
% MW_spv = 4; % 'kg/kmol'
% R_u = 8314.3; % 'J/kmol-K'
% T_T = 298; % 'K'
% T_crit_ox = 309.57; % 'K' % critical temperature of N2O
% V_T = 0.043; % 'm^3'
% 
% %% Derived Variables
% 
% A_inj = 0.25*pi*D_inj^2;
% n_oxt = m_oxt/MW_ox;
% n_spv = m_spv/MW_spv;
% n_oxl = (n_oxt*R_u*T_T - P_crit_ox(T_T)*V_T) / (-P_crit_ox(T_T)*V_mol_oxl + R_u*T_T); % initial N2O liquid [kmol]
% n_oxv = P_crit_ox(T_T)*(V_T - V_mol_oxl*n_oxl) / (-P_crit_ox(T_T)*V_mol_oxl + R_u*T_T); % initial N2O gas [kmol]

motor1 = motor('motor1');

%% State function initial condition and solution

state_0 = [motor1.n_oxv,motor1.n_oxl,motor1.T_T]; % kmol/s, kmol/s, K
tspan = 0:0.0005:5; % s

[t,state] = ode45(@(t,state) CV1_function(t,state,motor1),tspan,state_0);

%% State function definition
function state_dot = CV1_function(t,state,mi)  % mi = motor_instance
    
    P_C = 5; % bar

    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
    
    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez
    V_mol_oxl = mi.V_mol_oxl(T_T);
    P_crit_ox = mi.P_crit_ox(T_T);
    dP_crit_ox_dT = mi.dP_crit_ox_dT(T_T);
    CV_oxv = mi.CV_oxv(T_T);
    CV_spv = mi.CV_spv(T_T);

    C_d = mi.C_d;
    P_losses = mi.P_losses;
    N_inj = mi.N_inj;
    m_T = mi.m_T;
    MW_ox = mi.MW_ox;
    R_u = mi.R_u;
    T_T = mi.T_T;
    V_T = mi.V_T;
    A_inj = mi.A_inj;
    n_spv = mi.n_spv;
    
    P_T = (n_oxv + n_spv)*R_u*T_T/(V_T-n_oxl*V_mol_oxl);

    f1 = -C_d*N_inj*A_inj*sqrt((P_T-P_losses-P_C)/(2*MW_ox*V_mol_oxl));
    f2 = -V_mol_oxl*P_crit_ox;
    f3 = ((V_T - n_oxl*V_mol_oxl)*dP_crit_ox_dT - n_oxl*P_crit_ox*dV_mol_oxl_dT);
    f4 = R_u*n_oxv;
    f5 = R_u*T_T;
    f6 = m_T*c_P_T + n_oxl*CP_oxl + n_oxv*CV_oxv + n_spv*CV_spv;
    f7 = R_u*T_T - mi.deltaH_oxv(T_T);
    f8 = P_T*V_mol_oxl;

    M = [ 1 1 0; f5 -f2 f4-f3; f7 f8 -f6 ];
    V = [ f1; 0; 0 ];

    state_dot = linsolve(M,V);
end