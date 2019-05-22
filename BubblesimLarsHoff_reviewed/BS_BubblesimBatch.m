function BS_BubblesimBatch
% function BS_BubblesimBatch
%
% Simulate oscillation of bubble in acoustic field
%
% Calculate particle response by simulating a differential equation.
% Enter parameters manually

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

global BubblesimPath;
global nano micro milli centi kilo Mega;

ResultDirectory= sprintf('%s/Results', BubblesimPath );
cd(ResultDirectory);


%=== SIMULATION PARAMETERS =================================

%--- ODE solver ---
ODEsolvers=  {'ode45'
	      'ode15s'};

%--- Simulation models ---
models= {'Rayleigh'
	 'ModifiedRayleigh'
	 'Trilling'
	 'Keller' };

%--- Thermal damping options ---
thermaldampings= {'isothermal'
		  'adiabatic'
		  'pulsefrequency'
		  'resonancefrequency' };

%--- Pulse envelope ---
envelopes = {{'LoadFile',    0   }
             {'rectangular', 0   }
	     {'costapered',  0.2 }
	     {'triangular',  0   }
	     {'hanning',     0   }
	     {'hamming',     0   }
	     {'blackman',    0   }};


%--- Select model and solver ---
theory.ODEsolver     = ODEsolvers{2};
theory.model         = models{2};
theory.thermaldamping= thermaldampings{4};
theory

%--- Pulse parameters ---
envelope= envelopes{3};
envelope

invert= 0;
A = [ 100:100:800 1000 1200 1500]*kilo; % [Pa]  Pulse amplitude
Nc= [ 64 ];                             %       No. of pulse cycles
f0= [ 2 3 ]*Mega;               % [Hz]  Pulse center frequency
a0= [ 1.0:0.5:5.5 6:1:15 ]/2*micro;     % [m]   Particle radius

n=0;
for k2= 1:length(A );
  for k3= 1:length(Nc);
    for k4= 1:length(f0);
      for k5= 1:length(a0);
	n=n+1;
	parameter(n).invert  = invert;
	parameter(n).envelope= envelope;
	parameter(n).A = A(k2);
	parameter(n).Nc= Nc(k3);
	parameter(n).f0= f0(k4);
	parameter(n).a0= a0(k5);
      end
    end
  end
end
N= n;
drawnow
disp(sprintf('No. of parameter combinations: %d',N))
keyboard

%=== CALL CALCULATION FUNCTION =============================
for n=1:N
  [particle,pulse,simulation,graph]= CalculateBubblesim(parameter(n),theory);

  message= sprintf('Simulation %3d of %3d. Saved to %s', ...
                        n, N, graph.resultfile );
  WriteMessage(message);
  disp(message);
  drawnow;
end

cd ('..');
return


%=== FUNCTION: CALCULATE =======================================================
function [particle,pulse,simulation,graph]= CalculateBubblesim( parameter, theory );
global nano micro milli centi kilo Mega;

%==== PARAMETERS ===============================================================

%--- Particle ---
particle.a0 = parameter.a0;   % [m]   Bubble diameter
particle.ds = 4.0e-9;         % [m]   Shell thickness
particle.Gs = 50e6;           % [Pa]  Shell shear modulus
particle.es = 0.8;            % [Pas] Shell shear viscosity
particle.liquid.name='Water'; %       Surrounding liquid: 'Water' or 'Blood'

particle = PhysicalConstants(particle);  % Fixed physical constants

%--- Pulse ---
pulse.A  = parameter.A;
pulse.Nc = parameter.Nc;
pulse.f0 = parameter.f0;
pulse.fs = 100e6;
pulse.invert=parameter.invert;

pulse.source = '';
pulse.envelope.command= parameter.envelope;
pulse.envelope.name   = parameter.envelope{1};

%--- Simulation ---
simulation.solver.command        = theory.ODEsolver;
simulation.model.ode             = theory.model;
simulation.thermaldamping.command= theory.thermaldamping;
simulation.displayprogress       = 0;

%--- Plotting ---
graph.plotlinear = 1;
graph.include    = [1 1 1 0 1 0];
graph.fmax = min([ 4*pulse.f0, pulse.fs ]);
graph.tmax = 2*pulse.Nc/pulse.f0;


%=== DEFINE PULSE ==============================================================
[graph ]                 = DefineGraphs(graph);
[particle, pulse, graph] = DefinePulse(particle, pulse, graph);
[linear]                 = LinearOscillation (particle, pulse(1).t, pulse(1).p );

%=== SIMULATE RESULT ===========================================================
if (pulse(1).invert)
  N=2;
  simulation(2)=simulation(1);
else
  N=1;
end
 
for k=1:N
  tic;
  [t,a,ps] =SimulateOscillation (particle, pulse(k), simulation(k) );
  [tr,pr,fs]= ConstantSampleRate( t, ps, pulse(k).fs );

  simulation(k).t = t;    % [s]    Time vector from ODE solver,  uneven sampling
  simulation(k).a = a;    % [m]    Radius and velocity, uneven sampling
  simulation(k).p = ps;   % [Pa]   Scattered sound pressure, uneven sampling
  simulation(k).fs= fs;   % [1/s]  Sample rate, after resampling to constant rate
  simulation(k).tr= tr;   % [s]    Time vector, resampled to constant rate
  simulation(k).pr= pr;   % [Pa]   Scattered pressure, resampled to constant rate

  simulation(k).etime = toc;
end

%=== PLOT RESPONSE =============================================================
PlotInitialPulses ( particle, pulse, linear, graph );
PlotSimulation ( particle, pulse, linear, simulation, graph );
drawnow;

%=== SAVE RESULT ===============================================================
save ( graph.resultfile , 'particle', 'pulse', 'simulation' );


return














