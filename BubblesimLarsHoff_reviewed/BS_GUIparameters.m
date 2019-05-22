function [particle,pulse,simulation,graph] =BS_GUIparameters(particle,pulse,simulation,graph,action) 
% function [particle,pulse,simulation,graph] =BS_GUIparameters(particle,pulse,simulation,graph,action) 
%
% Bubblesim: Simulate bubble response
% Read or write particle and pulse parameters from/to GUI
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

global nano micro milli centi kilo Mega

BS_WriteFunctionname

%=== IDENTIFY HANDLES =============================================

%--- Particle
input.odesolver     = findobj( 'Tag','ODEsolver'         );   
input.ode           = findobj( 'Tag','ODE'               );   
input.radius        = findobj( 'Tag','Radius'            ); 
input.shellthickness= findobj( 'Tag','ShellThickness'    );
input.shellG        = findobj( 'Tag','ShellShearModulus' ); 
input.shelleta      = findobj( 'Tag','ShellViscosity'    ); 
input.liquid        = findobj( 'Tag','Liquid'           ); 

%--- Pulse
input.pulsetype= findobj('Tag','Pulsetype' );
input.A        = findobj('Tag','Amplitude' );
input.Nc       = findobj('Tag','Ncycles'   );
input.f0       = findobj('Tag','f0'        );
input.fs       = findobj('Tag','fs'        );

%--- Simulation
input.invert         = findobj('Tag','Inversion'      );
input.thermaldamping = findobj('Tag','ThermalDamping' );

%--- Plotting
input.plotlinear     = findobj('Tag','PlotLinear'     );
input.displayprogress= findobj('Tag','DisplayProgress');             
input.plot(1)        = findobj('Tag','PlotIncoming'   );
input.plot(2)        = findobj('Tag','PlotScattered'  );
input.plot(3)        = findobj('Tag','PlotRadius'     );
input.plot(4)        = findobj('Tag','PlotVelocity'   );
input.plot(5)        = findobj('Tag','PlotSpectra'    );
input.plot(6)        = findobj('Tag','PlotTransfer'   );

switch action
 %======== READ PARAMETERS FROM GUI ========================
 case 'read',  
  simulation.solver        = BS_GetUserdata( input.odesolver      );
  simulation.model         = BS_GetUserdata( input.ode            );
  simulation.thermaldamping= BS_GetUserdata( input.thermaldamping );   

  particle.a0= micro*str2num( get( input.radius,        'string') ); 
  particle.ds= nano *str2num( get( input.shellthickness,'string') );
  particle.Gs= Mega *str2num( get( input.shellG,        'string') ); 
  particle.es=       str2num( get( input.shelleta,      'string') ); 
  particle.liquid=  BS_GetUserdata(input.liquid);
  particle =        BS_PhysicalConstants(particle);
   
  pulse(1).envelope= BS_GetUserdata(input.pulsetype);
  pulse(1).A       = Mega*str2num( get( input.A , 'string' ) );
  pulse(1).Nc      =      str2num( get( input.Nc, 'string' ) );
  pulse(1).f0      = Mega*str2num( get( input.f0, 'string' ) );
  pulse(1).fs      = Mega*str2num( get( input.fs, 'string' ) );
  pulse(1).invert  = get( input.invert, 'value');
  
  simulation.displayprogress= get( input.displayprogress,  'value');  
  graph.plotlinear          = get( input.plotlinear,       'value');  
  for k=1:length(input.plot)
    graph.include(k) = get( input.plot(k), 'value');
  end
   
%=== WRITE PARAMETERS TO GUI ===============================
 case 'write'
  set (input.radius,         'string', particle.a0 *1e6  );
  set (input.shellthickness, 'string', particle.ds *1e9  );
  set (input.shellG,         'string', particle.Gs *1e-6 );
  set (input.shelleta,       'string', particle.es       );  
  set (input.A,              'string', pulse(1).A/Mega   );
  set (input.Nc,             'string', pulse(1).Nc       );
  set (input.f0,             'string', pulse(1).f0/Mega  );
  set (input.fs,             'string', pulse(1).fs/Mega  );
end

%== DEBUGGING: Display values ========================================
% $$$ input
% $$$ graph
% $$$ particle
% $$$ pulse
% $$$ simulation

return








