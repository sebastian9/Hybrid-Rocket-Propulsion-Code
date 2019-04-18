function [value, isterminal, direction] = stopIntegrationEvent(t, state)
    %% stops ode integration when the max height is reached
    if (state(2) <= 0)% && state(1) <= 0) % Linear momentum in z direction is zero
        value = 0; % when value = 0, an event is triggered
%     elseif (state(4) >= env.Motor.R_p) % Linear momentum in z direction is zero
%         value = 0; % when value = 0, an event is triggered
    else
        value =1;
    end
    isterminal = 1; % terminate after the first event
    direction = 0; % get all the zeros
end