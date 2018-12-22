function a = extract_data ( state,t)
% Function run the states through the model again to extract the needed
% internal parameters. This has to be done becasue ode45 does not allow for
% the extraction of internal parameters
global env
DTF_data =rocket(init_rocket());
motor_init( DTF_data );
env = environement(1400, 25, 86000, DTF_data );
%%
DTF_data.X = state(1,1:3)';
DTF_data.Q = state(1,4:7)';
DTF_data.P = state(1,8:10)';
DTF_data.L = state(1,11:13)';
DTF_data.time= t(1);
%%
ilast = length(t);
a = 0; %zeros(ilast,3);
%%
for i = 2:ilast
DTF_data.time= t(i);
DTF_data.deltat = t(i)-t(i-1);
DTF_data.X = state(i,1:3)';
DTF_data.Q = state(i,4:7)';
DTF_data.P = state(i,8:10)';
DTF_data.L = state(i,11:13)';
burn_data(DTF_data);
DTF_data.Xdot = (state(i,1:3)' - state(i-1,1:3)')/DTF_data.deltat;
DTF_data.Qdot = (state(i,4:7)' - state(i-1,4:7)')/DTF_data.deltat ;
DTF_data.Pdot = (state(i,8:10)' - state(i-1,8:10)')/DTF_data.deltat;
DTF_data.Ldot = (state(i,11:13)' -state(i-1,11:13)')/DTF_data.deltat;
Rmatrix= quat2rotm(DTF_data.Q');
RA = Rmatrix*env.RA0';
X = DTF_data.X;
Xdot = DTF_data.Xdot;
L = DTF_data.L;
if(norm(X) < DTF_data.Rail)
W = [0, 0, 0]';
else
W = env.W;
end
CnXcp = DTF_data.CnXcp;
Cn= CnXcp(1);
Xcp= CnXcp(2);
Cda = CnXcp(3); % Damping moment coefficient
zeta = CnXcp(4); % Damping ratio
Ssm = CnXcp(5); % Static stability margin
Ssm_B = CnXcp(6); % Static stability margin without body lift correction
Ccm = CnXcp(7); % Corrective moment coeff
invIbody = DTF_data.Ibody\eye(3); %inv(DTF.Ibody); inverting matrix
omega = Rmatrix*invIbody*Rmatrix'*L;
Vcm = Xdot + W;
Xstab = Xcp- DTF_data.Xcm;
omega_norm = normalize(omega); %normalized
Xprep =Xstab*sin(acos(dot(RA,omega_norm))); % Prependicular distance between omaga and RA
Vomega = Xprep *cross(RA,omega);
V = Vcm + Vomega; % approxamating the velocity of the cop
Vmag = norm(V);
Vnorm = normalize(V);
alpha = acos(dot(Vnorm,RA));
DTF_data.alpha = alpha;
%% Log Data
logData(DTF_data.alpha, DTF_data.Cd, Cda, DTF_data.Xcm, DTF_data.Mass, Vmag, Xcp, zeta, Ssm, Ssm_B, Ccm, t(i));
%logData(DTF_data.X(3),DTF_data.Cd,t(i)); % Eg DTF.Cd for drag norm(Xdot)/env.C
% log(i,1) = norm(DTF_data.Xdot);
% log(i,2) = DTF_data.Re;
% log(i,3) = t(i);
end
end
