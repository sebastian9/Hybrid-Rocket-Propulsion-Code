%% Paraffins
% All environment variables are defined here, this is 
% the master class passed to the solver
classdef Paraffins<handle
    properties
        MW
        deltaH
        T_m
        T_b
        T_v
        rho_s
        rho_l
        deltaH_F
        deltaH_V
        mu
        TC_l
        CP_s
        CV_l
    end
    methods
        function obj = Paraffins(paraffin_file)
            if nargin > 0
                
            end
        end
    end
end