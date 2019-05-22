function [particle,pulse,linear,simulation,graph]=BS_BubblesimCallback(action)
% function [particle,pulse,linear,simulation,graph]=BS_BubblesimCallback(action)
%
% Bubblesim: Simulate bubble response
% Callback functions
% Display and/or calculate bubble response
%
% action: 'display':   Display pulse
%         'calculate': Simulate result

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname
BS_WriteMessage('');
BS_WriteMessage('No warnings', 'warning');

switch(action)
 case 'display',
  [particle, pulse, linear, simulation, graph] = DisplayPulse;
 case 'calculate',
  [particle, pulse, linear, simulation, graph] = DisplayPulse;
  [particle, pulse, linear, simulation, graph] = ...
               BS_SimulateResult(particle, pulse, linear, simulation, graph);
 otherwise,
  error('CalculateBubblesim: Unknown calculation option')
end
return


%=== DISPLAY PULSES ========================================
function [particle,pulse,linear,simulation,graph] = DisplayPulse
BS_WriteFunctionname
BS_WriteMessage('Plotting input pulses')

particle  = [];  % Particle parameters
pulse     = [];  % Incoming pulse parameters
simulation= [];  % Simulation results
graph     = [];  % Plotting parameters

BS_WriteMessage('Reading input parameters')
[particle,pulse,simulation,graph]= ...
          BS_GUIparameters(particle,pulse,simulation,graph,'read');

BS_WriteMessage('Defining input pulse')
[particle, pulse, graph]= BS_DefinePulse(particle,pulse,graph);

graph.tmax = 0.5*max(pulse(1).t);   % Display settings. Can be zoomed in or out
graph.fmax = 4*pulse(1).f0;

graph = BS_DefineGraphs(graph);

[particle,pulse,simulation,graph]= ...
          BS_GUIparameters(particle,pulse,simulation,graph,'write');

BS_WriteMessage('Calculating linear results ')
[linear]= BS_LinearOscillation (particle, pulse(1).t, pulse(1).p );

BS_WriteMessage('Plotting input pulses')
BS_PlotInitialPulses ( particle, pulse, linear, graph );

set (graph.figure, 'Name', sprintf('Results  %s' , graph.title) );

BS_WriteMessage('Input pulses ready')
drawnow

return


%=== SIMULATE BUBBLE OSCILLATION ===========================
function [particle, pulse, linear, simulation, graph] = ...
   BS_SimulateResult(particle,pulse,linear,simulation,graph);

MessageWindow= gcf;
set (MessageWindow, 'Pointer', 'Watch'); % Display waiting indicator
drawnow;

%=== Simulate response =====================================
if (pulse(1).invert)
  N=2;
  simulation(2)=simulation(1);
else
  N=1;
end
 
for k=1:N
  tic
  BS_WriteMessage(sprintf('Simulating ODE numerically. Pulse %d', k) );
  [t, a, ps]= BS_SimulateOscillation (particle, pulse(k), simulation(k) );
  [tr,pr,fs]= BS_ConstantSampleRate( t, ps, pulse(k).fs );

  simulation(k).t = t;    % [s]    Time vector from ODE solver,  uneven sampling
  simulation(k).a = a;    % [m]    Radius and velocity, uneven sampling
  simulation(k).p = ps;   % [Pa]   Scattered sound pressure, uneven sampling
  simulation(k).fs= fs;   % [1/s]  Sample rate, after resampling to constant rate
  simulation(k).tr= tr;   % [s]    Time vector, resampled to constant rate
  simulation(k).pr= pr;   % [Pa]   Scattered pressure, resampled to constant rate

  simulation(k).etime = toc;
end

set (MessageWindow, 'Pointer', 'Arrow'); % Remove waiting indicator

%=== Plot results ==========================================
BS_WriteMessage('Plotting results')
BS_PlotSimulation ( particle, pulse, linear, simulation, graph );

save ( graph.resultfile , 'particle', 'pulse', 'simulation' );
BS_WriteMessage(sprintf('Finished. Results saved to file %s', graph.resultfile) );

return



