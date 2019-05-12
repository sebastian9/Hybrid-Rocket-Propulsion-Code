% The aerodynamic model is taken from the one described in 
% Chowdhury, S. M. (2012). DESIGN AND PERFORMANCE SIMULATION OF A HYBRID 
% SOUNDING ROCKET. Table 5-2 p.177

% Compressibility factor Eq 5-20 p. 160
if M >= 1
    beta = sqrt(M^2 -1);
else
    beta = sqrt(1- M^2);
end

% Supersonic Nose pressure drag (ogive & cone) Eq. 5-22 p.161
C_D_Nose_sup = 2.1*sin(epsilon_n)^2 + sin(epsilon_n)/(2*beta); % + dC_Ndalpha_nose * alpha_T^2
% Where:
% epsilon_n is the nose half angle
% beta is the compressibility factor
% alpha_T is the total angle of attack

% Subsonic Nose pressure drag (ogive & cone) Eq. 5-22 p.161
C_D_Nose_sub = C_D_Nose_sup * M ^ ((4/(1+k))*(1-sin(epsilon_n)/2)/sin(epsilon_n));
% Where:
% k is the specific heat ratio of air

% Cylindrical Fuselage Segment Pressure Drag
% Body subsonic segment body drag coefficient at angle of attack
% Eq-23 p.161 referred to Hopkins [58]
C_D_Body = 2*(k2-k1)*alpha_T^2*S_0/V_B^(2/3) + ...
    (2*alpha_T^3/V_B^(2/3)) * eta*r_body*C_D_c*(L_body - x_0);
% Where:
% (k2-k1) is the Munk apparent mass factor which is a function of the body fineness ratio.
% S_0 is the body cross sectional area
% V_B is the total body volume
% C_D_c is the steady-state cross-flow drag on an infinite cylinder
% x_0 is the body station where the flow ceases to be potential
% eta is the ratio of drag on a finite cylinder to the drag on an infinite cylinder
% r_body is the body radius
% L_body is the total body length

% Base drag (Boat tail) drag Eq-28
if M < 1
    C_D_Base = 0.12 + 0.13*M^2;
else
    C_D_Base = 0.12/M;
end
% Eq. 29
if motor_on
    C_D_Base = C_D_Base * (A_Base-A_noz)/(A_Base);
end

% Fin pressure drag Eq 5-27 for squared LE
if M < 1
	C_D_stag_cyl = .85*(1+.25*M^2+(M^4)/40);
else
	C_D_stag_cyl = .85*(1.84 - .76*M^-2 + .166*M^-4 + .035*M^-6);
end
C_D_Fin = N_Fin * (C_D_stag_cyl + C_D_Base) * S_0 / (b_Fin*t_fin*cos(lambda_Fin));
% Where:
% N_Fin is the fin count
% b_fin is the fin span
% lambda_fin is the fin sweep angle
% t_fin is the fin thickness


% Skin Friction Drag Eq 5-30
% Compresibility corrected skin friction
% drag coefficient
if Re <= 10^4
	C_f = 0.0148;
elseif Re > 10^4 && Re < 51*(R_s/L_body)^(-1.039)
	C_f = (1.5*ln(Re-5.6))^-2;
else
	C_f = 0.032*(R_s/L)^.2;
end
if M < 1
	C_f_c = C_f*(1-0.1*M^2);
elseif M >=1 && Re > 51*(R_s/L_body)^(-1.039)
	C_f_c = C_f*(1+.18*M^2)^-2;
else
	C_f_c = C_f*(1+0.15*M^2)^-.58;
end
% Drag coefficient due to skin friction
C_D_F = C_f_c * (S_0/S_WET);
% Where:
% Re is ASSUMED to be the reynolds number
% R_S is the surface roughnes height

