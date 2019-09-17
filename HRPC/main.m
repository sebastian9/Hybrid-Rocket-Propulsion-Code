clear
%% Environment Setup

% Initializes an environment with a given temperature and pressure, as well as
% a Motor for the analysis, the motor file contains all the relevant
% variables in standard units. Except for Pressure, Temperature and the ratio of specific
% heats, defined in the combustion chamber function.
env = Environments(100000,300,'motor1'); % Pa % K % motorfile ; % T = 293 Chelaru

%% State function initial conditions

n_oxv_0 = env.Motor.n_oxv; % kmol
n_oxl_0 = env.Motor.n_oxl; % kmol
T_T_0 = env.Motor.T_T; % K
R_p_0 = env.Motor.R_p; % m
P_C_0 = env.P; %Pa
T_C_0 = env.T; % K
rho_C_0 = P_C_0/(287.058*T_C_0); % kg/m^3
M_ox_C_0 = 0; % kg
M_f_C_0 = 0; % kg
M_a_C_0 = rho_C_0*(pi*R_p_0^2)*env.Motor.L_g; % kg
state_0 = [n_oxv_0, n_oxl_0, T_T_0, R_p_0, P_C_0, T_C_0, rho_C_0, M_ox_C_0, M_f_C_0, M_a_C_0, 0];
% state_0 = [n_oxv_0, n_oxl_0, T_T_0];

%% State Function Setup

% Initializes the ODE solver, it integrates the HRPC_function
% state equation through the given tspan, with the given initial conditions
% state_0, with the stop criteria defined in opts.
tspan = 0:0.005:50; % Default time span and step for the integration.
opts = odeset('Events', @stopIntegrationEvent); % Stop criteria for the ODE solver
[t,state] = ode45(@(t,state) HRPC_function(t,state,env),tspan,state_0,opts); % Initialization

%%  Plot results
PlotResults(t,state, env);
