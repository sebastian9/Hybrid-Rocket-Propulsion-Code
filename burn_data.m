function burn_data( DTF )
% Burn data function is executed at the start of each iteration in the
% solver. This updates the motor mass and inertias by assuming the
% impulse generated is proportional to the mass consumed
t= DTF.time;
T=DTF.T; % Thrust
tt=DTF.motordata(:,1);
TT=DTF.motordata(:,2);
ind= find(tt<t); % Finds index less then curretn time
ind=ind(end);
tt =[tt(1:ind); t]; % Burn uptill
TT= [TT(1:ind); T];
DTF.impulseGen = trapz(tt,TT); % impulse generated upto that point.
propM_used = DTF.propM_tot/DTF.Motor_impulse*DTF.impulseGen;
DTF.propM_current = DTF.propM_tot - propM_used; % Remaining prop mass
if(DTF.deltat == 0)
DTF.deltaMass =0 ;
else
% taking average for numerical stability
DTF.deltaMass = (DTF.propM_current - DTF.propM_prev)/ DTF.deltat; %1.2902/1.8;%
end
DTF.propM_prev = DTF.propM_current;
% New Inertias at time step
% IMP:The mass of the prop should already be updated in DTF to get correct Xcm
d = DTF.Xcm_prop - DTF.Xcm; % distance from Cm of propellent and???
prop_OD = DTF.prop_OD;
prop_h = DTF.prop_h;
Mass_prop = DTF.propM_current;
prop_density = DTF.prop_density;
prop_ID = sqrt(prop_OD^2 -4*Mass_prop/(pi*prop_h*prop_density ));
propIx = 0.5*Mass_prop*(prop_OD^2+prop_ID^2)/4;
propIy = Mass_prop/12*(3*(prop_ID^2+prop_ID^2)/4 + prop_h^2) + DTF.propM_current*(d);
DTF.Iprop = [propIx, 0, 0; 0, propIy, 0; 0, 0, propIy]; % propIy = propIz
end
