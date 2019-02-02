clear

env = Environments(101325,298,'motor1'); % Pa % K % motorfile

%% State function initial condition and solution

state_0 = [env.Motor.n_oxv,env.Motor.n_oxl,env.Motor.T_T]; % kmol/s, kmol/s, K
tspan = 0:0.0005:15; % s

opts    = odeset('Events', @stopEvent);
[t,state] = ode45(@(t,state) HRPC_function(t,state,env),tspan,state_0, opts);

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

%% Stop integration event
function [value, isterminal, direction] = stopEvent(t, state)
    value      = (state(2) <= 0)-0.5;
    isterminal = 1; % Stop the integration
    direction  = 0;
end