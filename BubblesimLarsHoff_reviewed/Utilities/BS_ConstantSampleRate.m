function [tr,xr,fs]= BS_ConstantSampleRate(t,x,fs)
% function [tr,xr,fs]= BS_ConstantSampleRate(t,x,fs)
%
% Resample to fixed sample rate fs 

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

fsi  = 1/min(diff(t));         %  Requested interpolation sample rate
fsmax= 2^15/(max(t)-min(t));   %  Max. allowed sample rate

fsi= max(fsi, fs   );          %  Interpolated rate, larger than fs
fsi= min(fsi, fsmax);          %  Interpolated rate, less than maximum

ti= [min(t):1/fsi:max(t)]';    %  Resample at constant rate
xi= interp1( t, x, ti, 'cubic'); 
  
%--- Anti-alias filter ---
SPT= exist('butter')==2;       %  Check whether Signal Processing TB is available
if SPT                         %  Filter only if Signal Processing TB
  fc= 0.4*fs;                  %  Cut-off frequency
  xf= BS_LPfilter(xi,fc,fsi);
else
  WriteMessage('Warning: No anti-alias filter. Signal Processing TB not found.', ...
	       'warning')
  xf= xi;                      % No low pass filtering available
end
tr= [min(ti):1/fs:max(ti)]';   % Sample down to fs. max(ti)<max(t)
xr= interp1( ti, xf, tr, 'cubic'); 

%--- Plot original and down-sampled results for comparison ---
%--- For debugging ---
% $$$  figure(3)
% $$$  clf
% $$$  plot( t, x,'b-' );
% $$$  hold on
% $$$  %plot( ti,xi,'g--');
% $$$  plot( tr,xr,'r:.');
% $$$  hold off

return


