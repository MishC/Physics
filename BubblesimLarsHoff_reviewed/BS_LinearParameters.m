function [delta, w0, K, kappa, eta]= BS_LinearParameters(particle, f)
% function [delta, w0]= BS_LinearParameters(particle, f)
%
% Calculate dimensionless damping constants
%  particle : Particle and environment parameters
%         f : Frequency [Hz]
%
% delta : Damping constants  delta= delta(radius, frequency, mechanism)
%     1 : Thermal conduction
%     2 : Shell viscosity
%     3 : Liquid viscosity
%     4 : Acoustic radiation
%     5 : Total damping
%    w0 : Angular resonance frequency [1/s]  w0=     w0(radius, frequency)
%     K : Bulk modulus [Pa]                   K=      K(radius, frequency) 
% kappa : Polytropic exponent             kappa=  kappa(radius, frequency) 
%   eta : Equivalent viscosities. See 'delta' for definitions.
% 
% Matrix organization
%   1: Particle radius
%   2: Frequency
%   3: Mechanism. See 'delta' for definition.


% Lars Hoff
% NTNU, Dept. of Telecommunications
% N-7491 Trondheim, Norway
% Revised 18 Sep 2000.

k0= find(f==0);

%--- Paticle and environment ---
a0 = particle.a0;
p0 = particle.p0;
rho= particle.rho;
c  = particle.c;
eL = particle.eL;

ds = particle.ds;
Gs = particle.Gs;
es = particle.es;

%--- Reshape to matrices ---
Nd= length(a0);
Nf= length(f);

w = repmat( 2*pi*f', Nd,  1);  
a = repmat( a0,       1, Nf);
ds= repmat( ds,       1, Nf);

delta = zeros(Nd,Nf,5);

%--- Thermal damping ---
[kappa, eTh, F, lD, X ]= ThermalDamping(particle,f );

%--- Resonance frequency ---
K  = kappa*p0 + 4*Gs*ds./a;   %  [Pa]   Particle bulk modulus
w02= 3*K./(rho*a.^2);         %  [1/s2] Angular resonance frequency
w0 = sqrt(w02);               %  [1/s]  Angular resonance frequency
W   = w./w0;                  %  [1]    Normalized frequency

%--- Radiation viscosity ---
ec = (rho*a.^3.*w.^2)./(4*c);

%--- Equivalent viscosities ---
eta(:,:,1)= eTh;                      % Thermal conduction
eta(:,:,2)= 3.*es.*ds./a;             % Shell viscosity
eta(:,:,3)= eL;                       % Liquid viscosity
eta(:,:,4)= (rho*a.^3.*w.^2)./(4*c);  % Acoustic radiation

%--- Dimensionless damping constants
Cd   = repmat(4./(w0.*rho.*a.^2), [1 1 4]); % Viscosity to damping constant
delta= Cd.*eta ;                            % Damping constant from viscosity

delta(:,:,5)= sum(delta(:,:,1:4), 3 );      % Total damping

return



%=== THERMAL DAMPING =================================================
function [kappa,eTh,F,lD,X]= ThermalDamping(particle,f);
%
% Calculate thermal damping constant and polythropic exponent
%	particle : Partticle and environment parameters
%	f [Hz]   : Frequency vector
%
%       F   [1]  : Thermal damping function
%       lD  [m]  : Thermal diffusion length
%       X        : a/lD

%--- Particle and environment parameters ---
a = particle.a0;    % [m]      Radius
p0= particle.p0;    % [Pa]     Ambient pressure
rg= particle.rg;    % [kg/m3]  Density of gas 
Cp= particle.Cp;    % [J/kg K] Heat capacity of gas at constant pressure
Kg= particle.Kg;    % [W/K m]  Thermal conductivity of gas
g = particle.gamma; % [1]      Cp/Cv

%--- Handle DC components separately ---
fz= (f==0);
kz= find(fz);
k = find(not(fz));  % Nonzero frequencies

%--- Define matrices ---
Nd= length(a);    %  Diameters
Nf= length(f);    %  Frequencies
Nk= length(k);    %  Nonzero frequencies
Nz= length(kz);   %  Zero frequencies

%--- Nonzero frequencies only ---
w = repmat( 2*pi*abs(f(k))', Nd,  1 );  % Nonzero, always positive
a = repmat( a,                1, Nk );

%--- Thermal damping function
lD = sqrt(Kg./(2*w*rg*Cp));      %  [m]  Thermal diffusion length
X  = a./lD;                      %  [1]  Radius/Diffusion length
Y  = (1+i)/2 *X;
F  = 1/g *(1+ 3*(g-1)./Y.^2 .*(Y.*coth(Y)-1) );

eTh  = (3*p0)./(4*w).*imag(1./F); % [Pas] Thermal damping viscosity
kappa= real(1./F);                % [1]   Polytropic exponent

%--- Low and high frequency limits ---
eTh0  = 1/20*(1-1/g)*(p0*rg*Cp)./(Kg)*a(:,1).^2; % Low frequency: Isothermal
eThInf= (3*p0)./(4*w)* 3*g*(g-1)./X;             % High frequency: Adiabatic

%--- Merge DC and AC components ---
       w(k)= w;                w(kz)=    0*ones(Nd,Nz);
 kappa(:,k)= kappa;      kappa(:,kz)=    1*ones(Nd,Nz);
   eTh(:,k)= eTh;          eTh(:,kz)= eTh0*ones(Nd,Nz);
     F(:,k)= F;              F(:,kz)=    1*ones(Nd,Nz);
      lD(k)= lD;              lD(kz)=  Inf*ones(Nd,Nz);
     X(:,k)= X;              X(:,kz)=    0*ones(Nd,Nz);
eThInf(:,k)= eThInf;    eThInf(:,kz)= Inf*ones(Nd,Nz);
    
return










