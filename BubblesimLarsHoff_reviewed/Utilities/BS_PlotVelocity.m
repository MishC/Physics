function ax= BS_PlotVelocity(Window,t,v,tmax,symbol )
% function ax= BS_PlotVelocity(Window,t,v,tmax,symbol )
%
%     Window : Figure window to plot in
%       t    : Time vector [s]
%       v    : Veocity vs. time [m/s]
%     symbol : Plot symbol code
%
%       ax   : Handle to plot
%
% Plot radial oscillation of particle 
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

%--- Find graph ---
if not (iscell(Window))
  ax=[];
  return
end

figure ( Window{1} );
subplot( Window{2} )
hold on

%--- Plot graph ---
plot(t, v, symbol )

xlabel ('Time [s]');
ylabel ('Velocity [m/s]')

%--- Set axes scales ---
if not(isempty(tmax))
  set (gca, 'xlim', [0 tmax ])
end

oldmax= axis; oldmax=oldmax(4);
newmax= max(BS_Scale125(abs(v)));
vmax  = max([oldmax; newmax]);
set (gca, 'ylim', vmax*[-1 1] )

zoom on

ax= gca;

return

