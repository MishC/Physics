function dxv= BS_Rayleigh( T, xv, flag, Q, parameter )
% function dxv= BS_Rayleigh( T, xv, flag, Q, parameter )
%
% Simulate Rayleigh-Plesset equation for gas encapsulated 
% in a shell
% Written for Matlab's ODE solvers
% Normalized radial displacement, pressure and time
%
% T     : Normalized time T = t*w0, w0= sqrt(p0/(rho*a0^2))
% x     : Radial displacement  a(t)= a0(1+x)
% flag  : Parameters for ODE solver, not used
% Q     : Normalized acoustic pressure: P= p(t)/p0
% parameter : Normalized visco-elastic parameters

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

qi= interp1( Q.t, Q.p, T, '*cubic' );   % Driving pressure
x = xv(1,:);  % Strain
dx= xv(2,:);  % Velocity

%--- Physical parameters ---
gs   = parameter.gs;
ns   = parameter.ns;
nL   = parameter.nL;
kappa= parameter.k;

%--- Pressure at bubble surface ---
[qL,q1,q2]= BS_SurfacePressure(x,dx,gs,ns,nL,kappa);

%--- ODE ---
q3 =  1+x;
ddx= -1./q3.*(3/2*dx.^2 +1+qi-qL );

%--- ODE as vector equation ---
dxv= [dx; ddx];

return



