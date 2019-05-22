function BS_PlotInitialPulses( particle, pulse, linear, graph );
% function BS_PlotInitialPulses( particle, pulse, linear, graph );
%
% Simulate oscillation of bubble in acoustic field
%
% Read parameters from menu window
% Define incoming acoustic pulse
% Define result windows
% Plot pulse, spectrum and linear responses
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname

%=== INCOMING PULSE AND SPECTRUM =====================================
if pulse(1).invert
  N=2;
else
  N=1;
end

for k=1:N
  BS_PlotPulse(graph.incoming, pulse(k).t, pulse(k).p, graph.tmax, pulse(k).A, graph.symbol.result{k});
end
BS_PlotSpectrum(graph.spectra, pulse(1).t, 1e-6*pulse(1).p, graph.fmax, pulse(1).fs, 'pulse', graph.symbol.incoming );


%=== LINEAR RESPONSE =================================================
if graph.plotlinear
   BS_PlotRadius   (graph.radius,    linear.t, linear.a(:,1), graph.tmax, particle.a0,          graph.symbol.linear );
   BS_PlotVelocity (graph.velocity,  linear.t, linear.a(:,2), graph.tmax,                       graph.symbol.linear );  
   BS_PlotPulse    (graph.scattered, linear.t, linear.p,      graph.tmax, [],                   graph.symbol.linear );
   BS_PlotSpectrum (graph.spectra,   linear.t, linear.p,      graph.fmax, pulse(1).fs, 'pulse', graph.symbol.linear );
end
BS_PlotSpectrum(graph.transfer, linear.f, linear.Hxp, graph.fmax, pulse(1).fs, 'spectrum', 'k-' );
BS_PlotSpectrum(graph.transfer, linear.f, linear.Hpp, graph.fmax, pulse(1).fs, 'spectrum', 'b-' );

return












