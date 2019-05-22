function [f,HdB]= BS_PlotSpectrum(Window,t,h,fmax,fs,Htype,symbol )
% function [f,HdB]= BS_PlotSpectrum(Window,t,h,fmax,fs,Htype,symbol )
%
%   Window : Handle to figure window and subplot
%        t : Time vector [Hz]
%        h : Impulse or frequency response [Rel. scale]
%     fmax : Max. frequency to plot 
%       fs : Sample rate used in resampling (Not actual rate)
%    Htype : Response type: 'pulse' or 'frequency'
%   symbol : Plot symbol code
%
%        f : Frequency vector [Hz]
%      AdB : Power amplitude [dB]
%
% Plot spectrum of pulse or impulse response 

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

%--- Find graph ---
if not (iscell(Window))
   f=0; HdB=0;
   return
end

%=== CALCULATE SPECTRUM ====================================
if strncmpi( Htype, 'pulse', 5 )   % Pulse or impulse response
  [ti,hi,fs]= BS_ConstantSampleRate(t, h, fs );
  Nfft= 2^nextpow2(2*length(ti));
  f  = [0:Nfft-1]'/Nfft *fs;
  Fh = fft(hi,Nfft)/fs;
else                               % Spectrum or transfer function
   Nfft= length(t);
   f= t; Fh= h; fs= 1;  
end

%--- Amplitude and phase spectrum ---
H  = abs(Fh).^2;       % Amplitude and phase 
phi= angle(Fh);
Np = 1:length(f)/2;    % Range containing information 
k  = find(H==0); H(k)=eps; 
HdB= 10*log10(H(Np));  % Amplitude in dB

%--- Plot amplitude ---
figure( Window{1} );
subplot( Window{2} )
hold on

plot( f(Np), HdB, symbol );     % Linear x-axis

xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]' )

%--- Axes maxima
if not(isempty(fmax))
  set(gca, 'xlim', [0 fmax])
end

if strcmp( get(gca,'ylimmode'),'manual')
   oldmax= axis;  oldmax= oldmax(4);
else
   oldmax=[];
end   
newmax= max([oldmax; ceil(max(HdB)/10)*10 ]);
set (gca, 'ylim', newmax-[60 0] );

grid off
zoom on

return


