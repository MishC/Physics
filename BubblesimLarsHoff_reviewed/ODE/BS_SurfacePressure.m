function [qL,q1,q2]= BS_SurfacePressure(x,dx,gs,ns,nL,kappa);
% function [qL,q1,q2]= BS_SurfacePressure(x,dx,gs,ns,nL,kappa);
%
% Pressure at the bubble surface
%
%   x : Radial strain
%  dx = dx/dt
%  gs : Normalized shell shear modulus
%  ns : Normalized shell shear viscosity
%  nL : Normalized liquid viscosity
%  kappa: Polytropic exponent
%
%  qL : Pressure at bubble surface
%  q1 = dqL/dx
%  q2 = dqL/ddx
%
%  Calcualate pressure at bubble surface
%  Boundary condition for ODE giving bubble motion
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

%--- Gas pressure ---
qg = (1+x).^(-3*kappa);
dqg= -3*kappa*(1+x).^(-3*kappa-1);

%--- Shell pressure ---
%--- Exponential shell model ---
x0= 1/8; eg= exp(-x/x0);  % Stiffness
x1= 1/4; en= exp(-x/x1);  % Viscosity

qs  = -12*( gs*x0*(1-eg) + ns   *en.*dx);
dqs1= -12*( gs*eg        - ns/x1*en.*dx);
dqs2= -12*ns*en;

%--- Pressure at bubble wall ---
qL =-4*nL*dx./(1+x)    +  qs + qg;
q1 = 4*nL*dx./(1+x).^2 +dqs1 +dqg;
q2 =-4*nL* 1./(1+x)    +dqs2;

return