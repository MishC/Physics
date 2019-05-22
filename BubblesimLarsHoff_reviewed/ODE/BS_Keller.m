function dxv= BS_Keller( T, xv, flag, pi, parameter )
% function dxv= BS_Keller( T, xv, flag, pi, parameter )
%
% Simulate Keller-Miksis equation for gas encapsulated in a shell
%
% Written for Matlab's ODE solvers
% Normalized radial displacement, pressure and time
%
% T     : Normalized time T = t*w0, w0= sqrt(p0/(rho*a0^2))
% x     : Radial displacement  a(t)= a0(1+x)
% flag  : Parameters for ODE solver, not used
% pi    : Normalized acoustic pressure: P= p(t)/p0
% parameter : Normalized visco-elastic parameters

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

qi= interp1( pi.t, pi.p, T, '*cubic' );   % Driving pressure
x = xv(1,:);  % Strain
dx= xv(2,:);  % Velocity

%--- Physical parameters ---
gs   = parameter.gs;
ns   = parameter.ns;
nL   = parameter.nL;
cn   = parameter.c;
kappa= parameter.k;
M  = dx/cn;

%--- Pressure at bubble surface ---
[qL,q1,q2]= BS_SurfacePressure(x,dx,gs,ns,nL,kappa);

%--- ODE ---
q3 = (1+x).*(1-M) - 1/cn*(1+x).*q2;
ddx= -1./q3.*(3/2*dx.^2.*(1-M/3) -(1+M).*(qL-1-qi)-1/cn*(1+x).*dx.*q1 );

%--- ODE as vector equation ---
dxv= [dx; ddx];

return














