function ax= BS_PlotSimulation( particle,pulse,linear,simulation,graph );
% function ax= BS_PlotSimulation( particle,pulse,linear,simulation,graph );
%
%      particle : Particle data
%         pulse : Incoming pulse
%        linear : Results of linear calculation
%    simulation : Simulation results
%         graph : Figure windows to plot in
%
%            ax : Handle to plot
%
% Bubblesim: Simulate bubble response
% Plot results of simulation
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname

for k=1:length(simulation)
   t = simulation(k).t;
   a = simulation(k).a(:,1);
   v = simulation(k).a(:,2);
   p = simulation(k).p;
   
   if not( isempty(a) )
      ax= BS_PlotRadius   (graph.radius,    t, a, graph.tmax, particle.a0, graph.symbol.result{k} );
      ax= BS_PlotVelocity (graph.velocity,  t, v, graph.tmax, graph.symbol.result{k}  );
   end
   if not(isempty(p))
      ax= BS_PlotPulse    (graph.scattered, t, p,  graph.tmax,  [],                  graph.symbol.result{k} );
      ax= BS_PlotSpectrum (graph.spectra,   t, p,  graph.fmax, pulse(k).fs, 'pulse', graph.symbol.result{k} );
   end
end

return
