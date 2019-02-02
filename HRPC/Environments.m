%% Environment
% All environment variables are defined here, this is 
% the master class passed to the solver
classdef Environments<handle
    properties
        R_u
        P
        T
        Motor
    end
    methods
        function obj = Environments(P,T,motor_file)
            if nargin > 0
                obj.P = P;
                obj.T = T;
                obj.R_u = 8314.3;
                obj.Motor = Motors(motor_file,obj);
            end
        end
    end
end