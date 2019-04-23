%% Motor
% Defines all the properties of the motor. The dynamics are continuously
% updated, any varibales need for functions should come from this class
classdef Motors <handle
    properties
        A_inj
        A_t
        A_e
        C_d_l
        c_P_T
        D_inj
        env
        L_g
        m_T
        m_oxt
        m_spv
        MW_ox
        MW_spv
        N_inj
        n_oxt
        n_oxl
        n_oxv
        n_spv
        P_losses
        P_C
        Q
        q
        R_p
        R_p_f
        rho_f
        T_crit_ox        
        T_T
        theta_noz
        V_T
    end
    methods
        % Constructor takes data from motor property file
        function obj = Motors(motor_file,env) % Motor file is the motor data
            if nargin > 0
                prop = readtable(motor_file,'Format','%s%f','HeaderLines', 1);
                prop = table2array(prop(1:end,2));
                obj.C_d_l = prop(1);
                obj.D_inj = prop(2);
                obj.m_oxt = prop(3);
                obj.m_spv = prop(4);
                obj.m_T = prop(5);
                obj.MW_ox = prop(6);
                obj.MW_spv = prop(7);
                obj.N_inj = prop(8);
                obj.P_losses = prop(9);
                obj.T_crit_ox = prop(10);
                obj.V_T = prop(11);
                obj.R_p = prop(12);
                obj.R_p_f = prop(13);
                obj.L_g = prop(14);
                obj.rho_f = prop(15);                
                obj.Q = prop(16); % Fuel energy density [J/kg]  
                obj.q = prop(17); % Heat loss through the chamber walls [W] 
                obj.A_t = prop(18);
                obj.A_e = prop(19);
                obj.theta_noz = prop(20);
                obj.P_C = env.P;
                obj.T_T = env.T;
                obj.env = env;
                % derived properties
                obj.A_inj = 0.25*pi*obj.D_inj^2;
                obj.n_oxt = obj.m_oxt/obj.MW_ox;
                obj.n_spv = obj.m_spv/obj.MW_spv;
                obj.n_oxl = (obj.n_oxt*obj.env.R_u*obj.T_T - obj.P_crit_ox(obj.T_T)*obj.V_T) ...
                    / (-obj.P_crit_ox(obj.T_T)*obj.V_mol_oxl(obj.T_T) + obj.env.R_u*obj.T_T); % initial N2O liquid [kmol]
                obj.n_oxv = obj.P_crit_ox(obj.T_T)*(obj.V_T - obj.V_mol_oxl(obj.T_T)*obj.n_oxt) ...
                    / (-obj.P_crit_ox(obj.T_T)*obj.V_mol_oxl(obj.T_T) + obj.env.R_u*obj.T_T); % initial N2O gas [kmol]
                obj.c_P_T = (4.8 + 0.00322*obj.T_T)*155.239; % 'J/(kg-K)' % specific heat of tank, Aluminium
            end
        end
        function [CV] = CV_oxv(obj,T)
            % heat capacity of N2O gas at constant pressure [J/(kmol*K)] coefficients, valid for Temp range [100 K - 1500 K]
            % D1 = 0.2934e5; D2 = 0.3236e5; D3 = 1.1238e3; D4 = 0.2177e5; D5 = 479.4;
            CV = 0.2934e5 + 0.3236e5*((1.1238e3/T)/sinh(1.1238e3/T))^2 + 0.2177e5*((479.4/T)/cosh(479.4/T))^2 - obj.env.R_u; 
        end
        function [CV] = CV_spv(obj,T)
            % heat capacity of He at constant pressure [J/(kmol*K)] coefficients, valid for Temp range [100 K - 1500 K]
            %C1 = 0.2079e5; C2 = 0; C3 = 0; C4 = 0; C5 = 0;
            CV = 0.2079e5 + 0*T + 0*T^2 + 0*T^3 + 0*T^4 - obj.env.R_u; %specific heat of He at constant volume [J/(kmol*K)]
        end
        function [H] = deltaH_oxv(obj,T)
            % heat of vaporization of N2O [J/kmol] coefficients, valid for Temp range [182.3 - 309.57 K]
            % J1 = 2.3215e7; J2 = 0.384; J3 = 0; J4 = 0;
            Tr = T/obj.T_crit_ox;
            H = 2.3215e7*(1 - Tr) ^ (0.384 + 0*Tr + 0*Tr^2);
        end
        function [T_C, CP_C, k_C] = NASACEA(obj, OF, P_C)
            fid = fopen ('test.inp', 'w');
            fprintf (fid, ' prob ro equilibrium \n\n');
            fprintf (fid, '  ! iac problem \n');
            fprintf (fid, ' o/f %f \n', OF);
            fprintf (fid, ' p,atm = %f \n', P_C/101325);
            fprintf (fid, ' pip %f \n', 1.00001*P_C/obj.env.P);
            fprintf (fid, ' reac \n');
            fprintf (fid, '   fuel  C32H66(a) wt%%=100 t,k=298.15 \n');
            fprintf (fid, '   oxid  N2O wt%%=100.  t,k=298.15 \n');
            fprintf (fid, ' output trace=1e-51 \n');
            fprintf (fid, ' end \n');
            fclose(fid);

            f = fopen('temp', 'w');
            fprintf (f, 'test\n');
            fclose (f);

            dos('FCEA2.exe < temp');

            text = fileread('test.out');
            text = extractAfter(text,strfind(text,'T, K')+16);
            text = extractBefore(text,8);
            T_C = str2double(text);
            text = fileread('test.out');
            text = extractAfter(text,strfind(text,'Cp, KJ/(KG)(K)')+17);
            text = extractBefore(text,7);
            CP_C = str2double(text);
            text = fileread('test.out');
            text = extractAfter(text,strfind(text,'GAMMAs')+17);
            text = extractBefore(text,7);
            k_C = str2double(text);
            
%             delete test.inp test.out temp
        end
    end
    methods (Static)
        function [CV] = CV_oxl(T) 
            % heat capacity of N2O liquid at constant pressure [J/(kmol*K) coefficients, valid for Temp range [182.3 K - 200 K]
            %E1 = 6.7556e4; E2 = 5.4373e1; E3 = 0; E4 = 0; E5 = 0;
            % specific heat of N2O liquid at constant volume, approx. same as at constant pressure [J/(kmol*K)]
            CV = 6.7556e4 + 5.4373e1*T + 0*T^2 + 0*T^3 + 0*T^4;
        end
        function [dPdT] = dP_crit_ox_dT(T)
            % vapor pressure of N2O [Pa] coefficients valid for Temp range [182.3 K - 309.57 K]
            % G1 = 96.512; G2 = -4045; G3 = -12.277; G4 = 2.886e-5; G5 = 2;
            % derivative of vapor pressure with respect T temperature
            dPdT = (--4045/(T^2) + -12.277/T + 2.886e-5*2*T^(2-1)) * exp(96.512 + -4045/T + -12.277*log(T) + 2.886e-5*T^2);
        end
        function [P] = P_crit_ox(T)
            % vapor pressure of N2O [Pa] coefficients valid for Temp range [182.3 K - 309.57 K]
            % G1 = 96.512; G2 = -4045; G3 = -12.277; G4 = 2.886e-5; G5 = 2;
            P = exp(96.512 + -4045/T + -12.277*log(T) + 2.886e-5*T^2);
        end
        function [V] = V_mol_oxl(T)
            % molar specific volume of liquid N2O [m^3/kmol] coefficients Q1 = 2.781; Q2 = 0.27244; Q3 = 309.57; Q4 = 0.2882;
            V = 0.27244.^(1+(1-T./309.57).^0.2882)./2.781;
        end
        function [dVdT] = dV_mol_oxl_dT(T)
            % molar specific volume of liquid N2O [m^3/kmol] coefficients Q1 = 2.781; Q2 = 0.27244; Q3 = 309.57; Q4 = 0.2882;
            dVdT = -(0.27244^((1 - T/309.57)^0.2882 + 1)*0.2882*log(0.27244)*(1 - T/309.57)^(0.2882 - 1))/(2.781*309.57);
        end
        function [R_p_dot] = R_p_dot(m_dot_ox,R_p,P_C)
%            Genevieve
%            a = 0.000155; n = 0.5; % Genevieve (2013) table 3.1
%            R_p_dot = a*(m_dot_ox/(pi*R_p^2))^n;
             % Chelaru, 2011 part 7
             a = 0.22e-4; n = 0.68; m = 0.07; l = 0.09;
             R_p_dot = a*(m_dot_ox/(pi*R_p^2))^n*P_C^m*(2*R_p)^l;
        end
    end
end
