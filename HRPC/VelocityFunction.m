function [result] = VelocityFunction(lambda, k_C, AtAe)
    % Chelaru, 2011, eq. 35
    a = (1-k_C)^(-1);
    b = (1-k_C)/(1+k_C);
    c = AtAe*((k_C+1)/2)^a;
    result = lambda-c*(1+b*lambda^2)^a;
end
