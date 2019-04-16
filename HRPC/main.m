clear

env = Environments(101325,286.5,'motor1'); % Pa % K % motorfile

%% State function initial condition and solution

n_oxv_0 = env.Motor.n_oxv; 
n_oxl_0 = env.Motor.n_oxl;
T_T_0 = env.Motor.T_T;
R_p_0 = env.Motor.R_p;
P_C_0 = env.P;
T_C_0 = env.T;
rho_C_0 = P_C_0/(287.058*T_C_0);
M_ox_C_0 = 0;
M_f_C_0 = 0;
M_a_C_0 = rho_C_0*(pi*R_p_0^2)*env.Motor.L_g;

state_0 = [n_oxv_0, n_oxl_0, T_T_0, R_p_0, P_C_0, T_C_0, rho_C_0, M_ox_C_0, M_f_C_0, M_a_C_0];
% state_0 = [env.Motor.n_oxv,env.Motor.n_oxl,env.Motor.T_T];
tspan = 0:0.05:15; % s

opts    = odeset('Events', @stopEvent);
[t,state] = ode45(@(t,state) HRPC_function(t,state,env),tspan,state_0,opts);

%%  Plot results 
figure(1), plot(t(:),state(:,3),'r','LineWidth',2),grid, ... 
    title('Temperature vs. Time'),... 
    xlabel('Time [s]'),... 
    ylabel('Temperature [K]');
figure(2), plot(t(:),state(:,1),'b',t(:),state(:,2),'g','LineWidth',2),grid, ... 
    title('kmol of N20 vs. Time'),... 
    xlabel('Time [s]'),... 
    ylabel('kmol of N2O [kmol]'),... 
    legend('kmol of N2O gas','kmol of N2O liquid');
figure(3), plot(t(:),state(:,4),'r','LineWidth',2),grid, ... 
    title('Port Radius vs. Time'),... 
    xlabel('Time [s]'),... 
    ylabel('Radius [m]');
figure(4), plot(t(:),state(:,5),'b','LineWidth',2),grid, ... 
    title('Chamber Pressure vs. Time'),... 
    xlabel('Time [s]'),... 
    ylabel('Pressure [Pa]');


%% Stop integration event
function [value, isterminal, direction] = stopEvent(t, state)
    %% stops ode integration when the max height is reached
    if (state(2) <= 0) % Linear momentum in z direction is zero
        value = 0; % when value = 0, an event is triggered
    else
        value =1;
    end
    isterminal = 1; % terminate after the first event
    direction = 0; % get all the zeros
end