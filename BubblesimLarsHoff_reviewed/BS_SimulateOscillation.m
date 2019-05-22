function [t, a, ps]= BS_SimulateOscillation( particle, pulse, simulation );
% function [t, a, ps]= BS_SimulateOscillation( particle, pulse, simulation );
%
% Bubblesim: Simulate bubble response
% Call ODE solvers
%
%  particle : Particle parameters
%  pulse    : Driving pulse parameters 
%
%  simulation: Result of simulation
%  t      : Time [s]
%  p      : Scattered pressure at bubble surface [Pa]
%  a(:,1) : Radius vs time [m]
%  a(:,2) : Velocity vs time [m/s]
% 

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname

%--- Estimate thermal damping parameters ---
switch simulation(1).thermaldamping.command
 case 'isothermal',  
  kappa = 1.0;
  eTh   = 0;
 case 'adiabatic'
  kappa = 1.4;
  eTh   = 0;
 case 'pulsefrequency'
  [delta, w0, K, kappa, eta]= ...
          BS_LinearParameters(particle, pulse(1).f0);  % Pulse center frequency
  eTh = eta(:,:,1);
 case 'resonancefrequency'
  [delta, w0, K, kappa, eta]= ...
	  BS_LinearParameters(particle, pulse(1).f0); 
  diff= inf;
  n=0;
  while diff>1e-6
    n=n+1;
    wold= w0;
    [delta, w0, K, kappa, eta]= ...
            BS_LinearParameters(particle, w0/(2*pi));   % Iterate to find resonance frequency
    diff = abs(w0-wold)/w0;
  end
  eTh = eta(:,:,1);
otherwise
 error('Unknown thermal damping model')
end;

%--- Driving pulse ---
t= pulse.t;
p= pulse.p;

%--- Paticle and environment ---
a0 = particle.a0;
p0 = particle.p0;
rho= particle.rho;
c  = particle.c;
eL = particle.eL + eTh;  % Effective viscosity increased from thermal conduction

ds = particle.ds;
Gs = particle.Gs;
es = particle.es;

%--- Normalized parameters ---
w0 = sqrt( p0/(rho*a0^2) ); % 1/s  Characteristic time scale

normalized.d = ds/a0;
normalized.nL= eL*w0/p0;
normalized.ns= es*w0/p0 *normalized.d;
normalized.gs=    Gs/p0 *normalized.d;
normalized.c = c/(a0*w0);
normalized.k = kappa;

q.t= t*w0;   % Normalized time
q.p= p/p0;   % Normalized pressure

%--- Simulate ODE ---
if simulation.displayprogress
  ODE.statistics= 'on';
  ODE.output    = 'odeplot';
  ODE.outputselect = [1];
else
  ODE.statistics= 'off';
  ODE.output    = [];
  ODE.outputselect = [];
end

options  = odeset('RelTol',     1e-6, ...
                  'AbsTol',    [1e-9 1e-9 ]', ...
		  'Stats',      ODE.statistics, ...
		  'OutputFcn',  ODE.output, ...
		  'OutputSel',  ODE.outputselect, ...
		  'Refine',     1, ...
		  'Vectorized','on' );
odefigure=figure;
[T,X]= feval( simulation.solver.command, ...
	      sprintf('BS_%s', simulation.model.ode), ...
	      [0 max(q.t)]', [0 0]', ...
          options, ...
	      q, normalized );
k = find(isnan(X)); X(k)= zeros(size(k));	% Remove NaNs from simulation
close(odefigure);

%--- Radiated pressure ---
qi.p= interp1( q.t, q.p, T, '*cubic' );   % Incoming pressure
qi.t= T;
qS= BubbleScatter( X, qi, normalized );

%--- Rescale to dimensional quantities ---
t     = T/w0;            % [s]    Time vector
a(:,1)= a0*(1+X(:,1));   % [m]    Radius vs. time
a(:,2)= a0*w0*X(:,2);    % [m/s]  Velocity vs. time
ps    = p0*a0*qS;

return



%=== RADIATED PRESSURE ===========================
function qS= BubbleScatter(X,Q,normalized);  
% 
% Calculate scattered pressure from oscillating bubble
%
% X  : Radial strain and velocity
% Q  : Incoming pressure
% normalized : Normalized visco-elastic parameters
%

t = Q.t;
q = Q.p;

x = X(:,1);
dx= X(:,2);

gs   = normalized.gs;
ns   = normalized.ns;
nL   = normalized.nL;
kappa= normalized.k;

%--- Pressure at bubble wall ---
[qL,q1,q2]= BS_SurfacePressure(x,dx,gs,ns,nL,kappa);

%--- Scattered pressure at the bubble wall ---
qS = (1+x).*(qL-1-q+1/2*dx.^2);

return









