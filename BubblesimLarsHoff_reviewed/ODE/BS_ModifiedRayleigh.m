function dxv= BS_ModifiedRayleigh( T, xv, flag, pi, parameter )
% function dxv= BS_ModifiedRayleigh( T, xv, flag, pi, parameter )
%
% Simulate modified Rayleigh-Plesset equation for gas encapsulated in a shell
% Equation modified by adding radiation damping term, dp/dt
%
% Written for Matlab's ODE solvers
% Normalized radial displacement, pressure and time
%
% T     : Normalized time T = t*w0, w0= sqrt(p0/(rho*a0^2))
% x     : Radial displacement  a(t)= a0(1+x)
% flag  : Parameters for ODE solver, not used
% pi    : Normalized acoustic pressure: P= p(t)/p0
% parameter : Normalized bubble and liquid parameters

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

%--- Pressure at bubble surface ---
[qL,q1,q2]= BS_SurfacePressure(x,dx,gs,ns,nL,kappa);

%--- ODE ---
q3 = (1+x) -1/cn.*(1+x).*q2;
ddx=-1./q3.*(3/2*dx.^2 -1/cn*(1+x).*q1.*dx +1+qi-qL );

%--- ODE as vector equation ---
dxv= [dx; ddx];

return








