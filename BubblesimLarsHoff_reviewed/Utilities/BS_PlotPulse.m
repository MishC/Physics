function ax= BS_PlotPulse(Window,t,p,tmax,A,symbol )
% function ax= BS_PlotPulse(Window,t,p,A,symbol )
%
%  Window : Handle to figure window and subplot
%       t : Time vector [s]
%       p : Pressure vs. time [Pa]
%    tmax : Max time to plot
%       A : Amplitude, scale maximum [Pa]
%  symbol : Plot symbol code
%
%      ax : Handle to plot
%
%
% Plot ultrasound pulse
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

%--- Find graph ---
if not (iscell(Window))
   ax=[];
   return
end
figure( Window{1} );
subplot( Window{2} )
hold on

%--- Plot pulse
plot(t, p, symbol)

xlabel ('Time [s]');
ylabel ('Pressure [Pa]')


%--- Set axis maxima
ax= axis;
if not(isempty(tmax))
  set(gca, 'xlim', [0 tmax])
end

if not(isempty(A))
  pm=A;           % Use specified y-axis limit
else
  if strcmp( get(gca,'ylimmode'),'manual')
    oldmax= ax(4);  
  else
    oldmax=[];
  end  
  pm= max([oldmax; BS_Scale125(abs(p)) ]);   % Find y-axis maximum
end   
set ( gca, 'ylim',  pm*[-1 1] )

ax= gca;

zoom on

return








