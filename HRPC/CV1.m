clear

motor1 = motor('motor1');

%% State function initial condition and solution

state_0 = [motor1.n_oxv,motor1.n_oxl,motor1.T_T]; % kmol/s, kmol/s, K
tspan = 0:0.0005:15; % s

opts    = odeset('Events', @stopEvent);
[t,state] = ode45(@(t,state) CV1_function(t,state,motor1),tspan,state_0, opts);

% Plot results 
figure(1), plot(t(:),state(:,3),'r','LineWidth',2),grid, ... 
    title('Temperature vs. Time'),... 
    xlabel('Time [s]'),... 
    ylabel('Temperature [K]');
figure(2), plot(t(:),state(:,1),'b',t(:),state(:,2),'g','LineWidth',2),grid, ... 
    title('kmol of N20 vs. Time'),... 
    xlabel('Time [s]'),... 
    ylabel('kmol of N2O [kmol]'),... 
    legend('kmol of N2O gas','kmol of N2O liquid');

%% State function definition
function state_dot = CV1_function(t,state,mi)  % mi = motor_instance
    
    % Curve fitted combustion chamber pressure [Pa]: 
    P_C = -2924.42*t^6 + 46778.07*t^5 - 285170.63*t^4 + 813545.02*t^3 - ... 
        1050701.53*t^2 + 400465.85*t + 1175466.2;

    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
    
    CP_oxl = mi.CV_oxl(T_T); % approximation from Fernandez
    V_mol_oxl = mi.V_mol_oxl(T_T);
    P_crit_ox = mi.P_crit_ox(T_T);
    dP_crit_ox_dT = mi.dP_crit_ox_dT(T_T);
    dV_mol_oxl_dT = mi.dV_mol_oxl_dT(T_T);
    CV_oxv = mi.CV_oxv(T_T);
    CV_spv = mi.CV_spv(T_T);

    C_d = mi.C_d_l;
    c_P_T = mi.c_P_T;
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

function [value, isterminal, direction] = stopEvent(t, state)
    value      = (state(2) <= 0)-0.5;
    isterminal = 1;   % Stop the integration
    direction  = 0;
end