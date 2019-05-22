function [particle, pulse, graph]= BS_DefinePulse(particle, pulse, graph);
% function [particle, pulse, graph]= BS_DefinePulse(particle, pulse, graph);
%
% particle : Particle parameters 
%    pulse : Incoming pulse data
%    graph : Graph and plotting parameters
%
% Simulate oscillation of bubble in acoustic field
% Define driving acoustic pulse
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

global nano micro milli centi kilo Mega

BS_WriteFunctionname

%=== DEFINE PULSE ===============================================================
if strncmpi( pulse.envelope.command{1}, 'load', 4 );
   pulse= LoadPulse(pulse);
else
   [pulse.t pulse.p]= BS_MakePulse( pulse.A, pulse.Nc, pulse.f0, pulse.fs, pulse.envelope.command );
   pulse.source='Synthetic pulse';
   pulse(2)    = pulse(1);
   pulse(2).p  =-pulse(1).p;  % Inverted incoming pulse
end

graph.title = sprintf('%s, %4.2f MHz, %5.3f MPa', ...
     pulse(1).source, pulse(1).f0/Mega, pulse(1).A/Mega ); 

graph.resultfile = BS_MakeFilename;

return


%=== LOAD PULSE FROM FILE ======================================
function [pulse] = LoadPulse(pulsedata)

[filename, pathname] = uigetfile('*.pls', 'Select pulse file');
load ( [pathname, filename], '-mat' );
cd (pathname)

%=== Calculate new parameter values
pulse.name   = strtok(filename, '.');
pulse.source = filename;   
pulse.envelope.command= 'Loaded';   

pulse.invert = 0;
pulse.Nc = 0;
pulse.f0 = pulse.f;

pulse.f = [];

return



