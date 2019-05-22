function linear= BS_LinearOscillation(particle,t,p);
% function linear= BS_LinearOscillation(particle,t,p);
%
% particle : Particle parameters
%       t  : time vector [s]
%       p  : Pressure pulse [Pa]
%
% linear : Results of linear model calculation
% Fields
%     t: Time [s]
%     x: Radial displacement [m]
%     a: Radius (1) [m] and velocity (2) [m/s]
%     p: Pressure [Pa]
%     f: Frequency [Hz]
%   Hxp: Transfer function pressure in - radius
%   Hpp: Transfer function pressure in - pressure out
%   hxp: Impulse response pressure in - radius
%   hpp: Impulse response pressure in - pressure out
%
% Simulate oscillation of bubble in acoustic field
% Calculate radial oscillation
% Linear model, frequency domain transfer functions
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname

%--- Physical parameters ---
a0 = particle.a0;
rho= particle.rho;

%--- Driving pulse ---
Ni = max(size(t));
reshape(t,Ni,1);
reshape(p,Ni,1);

%--- Pulse parameters
dt = t(2)-t(1);
fs = 1/dt;		% Sampling frequency

%--- FFT parameters
N = 2^nextpow2(2*Ni);
fb= fs/2* [1:N/2-1]'/(N/2);
f = [0;fb;fs/2;-flipud(fb)];  
w = 2*pi*f;
k0= find(f==0);

%--- Pulse spectrum
Fp = fft(p, N); 

%--- Transfer functions
[delta, w0, K, kappa, eta]= BS_LinearParameters(particle, f);
delta= delta(:,:,5);  % Total damping 
delta= delta';        % Reshape to column vectors. Ref.: 'LinearParameters'
w0   = w0';

%---Pressure to radius ----
W  = w./w0; 
Cx = 1./(rho*a0*w0.^2); 
Hxp= 1./( W.^2-1-i.*W.*delta );  Hxp(k0)=-1;
Fx = Cx.*Hxp.*Fp;

%--- Radiated pressure at unit distance ---
Hpp=-W.^2.*Hxp;     % [1]  Transfer function: Pressure in -> Pressure at surface
aN = 1;             % [m]  Unit distance
Fpl= a0/aN*Hpp.*Fp; %      Radiated pressure at unit distance     

%--- Time traces by inverse FT 
tl = dt*[0:1:N-1]';      % Time vector including padded zeros
xl = real(ifft(Fx));     % Radial oscillation
pl = real(ifft(Fpl));    % Radiated pressure
hxp= real(ifft(Hxp));    % Impulse response, radial amplitude
hpp= real(ifft(Hpp));    % Impulse response, radiated pressure

nz = [ N*(1-1/16):N];    % End of impulse response: Set to zero to remove oscillations
hpp(nz)= zeros(size(nz)); 

u=  [diff(xl)./diff(tl); 0];  % Radial velocity
al= [xl+a0 u];	              % Instantanous radius and velocity 

%--- Put results in structure
linear.t = tl;
linear.x = xl;
linear.a = al;
linear.p = pl;
linear.f = f;
linear.Hxp= Hxp;
linear.Hpp= Hpp;
linear.hxp= hxp;
linear.hpp= hpp;

return
