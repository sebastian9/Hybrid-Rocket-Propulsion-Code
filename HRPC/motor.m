%% Motor
% Defines all the properties of the motor. The dynamics are continuously
% updated, any varibales need for functions should come from this class
classdef motor <handle
    properties
        C_d_l
        D_inj
        P_losses
        N_inj
        m_T
        m_oxt
        m_spv
        MW_ox
        MW_spv
        R_u
        T_T
        T_crit_ox
        V_T
        A_inj
        n_oxt
        n_spv
        n_oxl
        n_oxv
        c_P_T
    end
    methods
        % Constructor takes data from motor property file
        function obj = motor(motor_file) % Motor file is the motor data
            if nargin > 0
                prop = readtable(motor_file,'Format','%s%f');
                prop = table2array(prop(1:end,2));
                obj.C_d_l = prop(1);
                obj.D_inj = prop(2);
                obj.P_losses = prop(3);
                obj.N_inj = prop(4);
                obj.m_T = prop(5);
                obj.m_oxt = prop(6);
                obj.m_spv = prop(7);
                obj.MW_ox = prop(8);
                obj.MW_spv = prop(9);
                obj.R_u = prop(10);
                obj.T_T = prop(11);
                obj.T_crit_ox = prop(12);
                obj.V_T = prop(13);
                % derived properties
                obj.A_inj = 0.25*pi*obj.D_inj^2;
                obj.n_oxt = obj.m_oxt/obj.MW_ox;
                obj.n_spv = obj.m_spv/obj.MW_spv;
                obj.n_oxl = (obj.n_oxt*obj.R_u*obj.T_T - obj.P_crit_ox(obj.T_T)*obj.V_T) / (-obj.P_crit_ox(obj.T_T)*V_mol_oxl(obj.T_T) + obj.R_u*obj.T_T); % initial N2O liquid [kmol]
                obj.n_oxv = obj.P_crit_ox(obj.T_T)*(obj.V_T - V_mol_oxl(obj.T_T)*obj.n_oxl) / (-obj.P_crit_ox(obj.T_T)*V_mol_oxl(obj.T_T) + obj.R_u*obj.T_T); % initial N2O gas [kmol]
                obj.c_P_T = (4.8 + 0.00322*obj.T_T)*155.239; % 'J/(kg-K)' % specific heat of tank, Aluminium
            end
        end
    end
    methods (Static)
        function [CV] = CV_oxl(T) 
            % heat capacity of N2O liquid at constant pressure [J/(kmol*K) coefficients, valid for Temp range [182.3 K - 200 K]
            E1 = 6.7556e4; E2 = 5.4373e1; E3 = 0; E4 = 0; E5 = 0;
            % specific heat of N2O liquid at constant volume, approx. same as at constant pressure [J/(kmol*K)]
            CV = E1 + E2*T + E3*T^2 + E4*T^3 + E5*T^4;
        end
        function [CV] = CV_oxv(T)
            % heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients
            % valid for Temp range [100 K - 1500 K]
            D1 = 0.2934e5; D2 = 0.3236e5; D3 = 1.1238e3; D4 = 0.2177e5; D5 = 479.4;
            %specific heat of N2O gas at constant volume [J/(kmol*K)]
            CV = D1 + D2*((D3/T)/sinh(D3/T))^2 + D4*((D5/T)/cosh(D5/T))^2 - R; 
        end
        function [CV] = CV_spv(T)
            % heat capacity of He at constant pressure [J/(kmol*K)] coefficients, valid for Temp range [100 K - 1500 K]
            C1 = 0.2079e5; C2 = 0; C3 = 0; C4 = 0; C5 = 0;
            CV = C1 + C2*T + C3*T^2 + C4*T^3 + C5*T^4 - R; %specific heat of He at constant volume [J/(kmol*K)]
        end
        function [H] = deltaH_oxv(T)
            obj.T_crit_ox = 309.57; % 'K'
            % heat of vaporization of N2O [J/kmol] coefficients, valid for Temp range [182.3 - 309.57 K]
            J1 = 2.3215e7; J2 = 0.384; J3 = 0; J4 = 0;
            Tr = T/obj.T_crit_ox;
            H = J1*(1 - Tr) ^ (J2 + J3*Tr + J4*Tr^2);
        end
        function [dPdT] = dP_crit_ox_dT(T)
            % vapor pressure of N2O [Pa] coefficients 
            % valid for Temp range [182.3 K - 309.57 K]
            G1 = 96.512; G2 = -4045; G3 = -12.277; G4 = 2.886e-5; G5 = 2;
            % derivative of vapor pressure with respect T temperature
            dPdT = (-G2/(T^2) + G3/T + G4*G5*T^(G5-1)) * exp(G1 + G2/T + G3*log(T) + G4*T^G5);
        end
        function [P] = P_crit_ox(T)
            % vapor pressure of N2O [Pa] coefficients 
            % valid for Temp range [182.3 K - 309.57 K]
            G1 = 96.512; G2 = -4045; G3 = -12.277; G4 = 2.886e-5; G5 = 2;
            P = exp(G1 + G2/T + G3*log(T) + G4*T^G5);
        end
        function [V] = V_mol_oxl(T)
            % molar specific volume of liquid N2O [m^3/kmol] coefficients
            Q1 = 2.781; Q2 = 0.27244; Q3 = 309.57; Q4 = 0.2882;
            V = Q2^(1+(1-T/Q3)^Q4)/Q1;
        end
    end
end
