function particle = BS_PhysicalConstants(particle)
% function particle = BS_PhysicalConstants(particle)
%
% particle : Particle and environment 
%
% Simulate oscillation of bubble in acoustic field
% Define values of physical constants
% Gas and liquid thermal and mechanical properties
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

%--- Fixed physical parameters ---
switch lower(particle.liquid.name)
case 'water',
 particle.p0 = 1.013e5;   % [Pa]     Ambient pressure
 particle.rho= 1000;      % [kg/m3]  Density of liquid
 particle.eL = 1.0e-3;    % [Pas]    Viscosity in liquid
 particle.c  = 1500;      % [m/s]    Speed of sound in liquid
case 'blood', % Ref. Angelsen: Waves, Signals, ... Table 2-1  
 particle.p0 = 1.013e5;   % [Pa]     Ambient pressure
 particle.rho= 1025;      % [kg/m3]  Density of liquid
 particle.eL = 4e-3;      % [Pas]    Viscosity in liquid
 particle.c  = 1570;      % [m/s]    Speed of sound in liquid
otherwise
 error ('Unknown liquid')
end

%--- Gas parameters. Values for air ---
particle.Kg   = 26.2e-3;   % [W/(mK)]   Thermal conductivity
particle.rg   = 1.161;     % [kg/m3]    Density of gas
particle.Cp   = 1.007e3;   % [J/(kgK)]  Heat capacity, const. pressure
particle.gamma= 1.4;       % [1]        Cp/Cv

return