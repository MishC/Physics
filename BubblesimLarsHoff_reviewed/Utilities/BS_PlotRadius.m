function ax= BS_PlotRadius(Window,t,a,tmax,a0,symbol )
% function ax= BS_PlotRadius(Window,t,a,a0,symbol )
%
%     Window : Figure window to plot in
%       t    : Time vector [s]
%       a    : Radius vs. time [m]
%       a0   : Equilibrium radius [m]
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
plot( t, a, symbol )

xlabel ('Time [s]');
ylabel ('Radius [m]')

%--- Set axes scales ---
if not(isempty(tmax))
  set (gca, 'xlim', [0 tmax ])
end

amax= ceil(max(a*1e6))*1e-6;
set (gca, 'ylim', [0 amax ])
zoom on;

ax= gca;

return
