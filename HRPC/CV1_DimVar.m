% Consideradas: A_inj C_d c_P_T CP_oxl CV_oxv CV_spv m_T MW_ox N_inj
% dP_crit_ox_dT P_crit_ox P_losses R_u T_T deltaH_oxv
% No consideradas: n_oxv n_oxl n_spv P_T P_CV_T V_mol_oxl dV_mol_oxl_dT

%% Primal Variables

D_inj = DimVar(0.002,'m');
C_d_l = 0.8;
% C_d_g = 0.5;
C_d = C_d_l; 
P_losses = DimVar(5,'bar');
R_u = DimVar(8314.3, 'J/kmol-K');
N_inj = 17;
m_T = DimVar(6.4882, 'kg');
MW_ox = DimVar(44.013, 'kg/kmol');

%% Derived Variables

A_inj = 0.25*pi*D_inj^2;

%% Perry's Chemical Engineers' Handbook Property Equations

% vapor pressure of N2O [Pa] coefficients % valid for Temp range [182.3 K - 309.57 K]
G1 = 96.512; G2 = -4045; G3 = -12.277; G4 = 2.886e-5; G5 = 2;
% critical temperature of N2O [K] 
T_crit_ox = DimVar(309.57,'K');
% heat of vaporization of N2O [J/kmol] coefficients, valid for Temp range [182.3 - 309.57 K]
J1 = 2.3215e7; J2 = 0.384; J3 = 0; J4 = 0;
% heat capacity of He at constant pressure [J/(kmol*K)] coefficients, valid for Temp range [100 K - 1500 K]
C1 = 0.2079e5; C2 = 0; C3 = 0; C4 = 0; C5 = 0;
% heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients, valid for Temp range [100 K - 1500 K]
D1 = 0.2934e5; D2 = 0.3236e5; D3 = 1.1238e3; D4 = 0.2177e5; D5 = 479.4;
% heat capacity of N2O liquid at constant pressure [J/(kmol*K) coefficients, valid for Temp range [182.3 K - 200 K]
E1 = 6.7556e4; E2 = 5.4373e1; E3 = 0; E4 = 0; E5 = 0;
% molar specific volume of liquid N2O [m^3/kmol] coefficients
Q1 = 2.781; Q2 = 0.27244; Q3 = 309.57; Q4 = 0.2882;

%% State function initial condition and solution

state_0 = [DimVar(0,'kmol/s'); DimVar(0,'kmol/s'); DimVar(298,'K')];
tspan = DimVar(0:0.0005:5,'s');

[t,state] = ode45(@CV1_function,tspan,state_0);

%% State function definition
function state_dot = CV1_function(t,state)
    
    P_C = DimVar(5,'bar');

    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
    
    n_to = m_loaded/MW2; % initial total N2O in tank [kmol]  
    V_mol_oxl = DimVar(Q2^(1+(1-Ti/Q3)^Q4)/Q1,'m^3/kmol'); % molar volume of liquid N2O
    P_crit_ox = DimVar(exp(G1 + G2/T_T + G3*log(T_T) + G4*T_T^G5),'Pa'); % initial vapor pressure of N20

    n_ox_v = P_crit_ox*(V - V_mol_oxl*n_to) / (-P_crit_ox*V_mol_oxl + R*T_T); % initial N2O gas [kmol] 
    n_ox_g = (n_to*R*T_T - P_crit_ox*V) / (-P_crit_ox*V_mol_oxl + R*T_T); % initial N2O liquid [kmol]
    
    c_P_T = DimVar((4.8 + 0.00322*T_T)*155.239,'J/(kg-K)'); % specific heat of tank, Aluminium
    
    CV_spv = C1 + C2*To + C3*To^2 + C4*To^3 + C5*To^4 - R; %specific heat of He at constant volume [J/(kmol*K)]
    CV_oxv = D1 + D2*((D3/To)/sinh(D3/To))^2 + D4*((D5/To)/cosh(D5/To))^2 - R; %specific heat of N2O gas at constant volume [J/(kmol*K)]
    % specific heat of N2O liquid at constant volume, approx. same as at constant pressure [J/(kmol*K)] % reduced temperature
    % heat of vaporization of N2O [J/kmol] % vapor pressure of N20 [Pa]
    dP_crit_ox_dT = (-G2/(To^2) + G3/To + G4*G5*To^(G5-1)) * exp(G1 + G2/To + G3*log(To) + G4*To^G5); % derivative of vapor pressure with respect to temperature
   
    % Given functions of temperature: Vhat_l = Q2^(1+(1-To/Q3)^Q4)/Q1; %molar specific volume of liquid N2O [m^3/kmol]
    CV_oxl = E1 + E2*To + E3*To^2 + E4*To^3 + E5*To^4;
    CP_oxl = CV_oxl; % approximation from Fernandez
    Tr = T_T/T_crit_ox;
    deltaH_oxv = J1*(1 - Tr) ^ (J2 + J3*Tr + J4*Tr^2);

    f1 = -C_d*N_inj*A_inj*sqrt((P_T-P_losses-P_C)/(2*MW_ox*V_mol_oxl));
    f2 = -V_mol_oxl*P_crit_ox;
    f3 = ((V_T - n_oxl*V_mol_oxl)*dP_crit_ox_dT - n_oxl*P_crit_ox*dV_mol_oxl_dT);
    f4 = R_u*n_oxv;
    f5 = R_u*T_T;
    f6 = m_T*c_P_T + n_oxl*CP_oxl + n_oxv*CV_oxv + n_spv*CV_spv;
    f7 = R_u*T_T - deltaH_oxv;
    f8 = P_T*V_mol_oxl;

    M = [ 1 1 0; f5 -f2 f4-f3; f7 f8 -f6 ];
    V = [ f1; 0; 0 ];

    state_dot = linsolve(M,V);

end

