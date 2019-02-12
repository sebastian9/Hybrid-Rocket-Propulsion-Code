syms Mxdot Mydot

Mx = 0;
mz = 0.1;
My = 0;
mx = 0.5;
my = 0.1;

[Mxdot, Mydot] = solve( [ Mxdot == mx - mz / ( 1 + (My+Mydot)/(Mx+Mxdot) ) ,...
    Mydot == my - mz / ( 1+(Mx+Mxdot)/(My+Mydot) ) ] , [Mxdot,Mydot] )

OF_m = mx/my