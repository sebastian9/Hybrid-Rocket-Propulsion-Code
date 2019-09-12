function state_dot = HRPC_function(t,state,env)
    %% State Variables
    n_oxv = state(1);
    n_oxl = state(2);
    T_T = state(3);
    R_p = state(4);
    P_C = state(5);
    T_C = state(6);
    rho_C = state(7);
    M_ox_C = state(8);
    M_f_C = state(9);
    M_a_C = state(10);
%     Pressure regressions for testing Fernadez code. Fernandez, 2009 Appendix H
%     P_C = -2924.42*t^6 + 46778.07*t^5 - 285170.63*t^ 4 + 813545.02*t^3 - ...
%         1050701.53*t^2 + 400465.85*t + 1175466.2; % Fernandez Motor Test 1
%     P_C = 95.92*t^6 - 2346.64*t^5 + 21128.78*t^4 - 87282.73*t^3 + ...
%         186675.17*t^2 - 335818.91*t + 3029190.03;  % Fernandez Motor Test 2
    %% Phase-Space equations
    state_dot(1:3,1) = OxidizerTank(n_oxv, n_oxl, T_T, P_C, env.Motor);
    state_dot(4:11,1) = CombustionChamber(state_dot(1), state_dot(2), R_p, T_C, P_C, rho_C, M_ox_C, M_f_C, M_a_C, env.Motor);
end
