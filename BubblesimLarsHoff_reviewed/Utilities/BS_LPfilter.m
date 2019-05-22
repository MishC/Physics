function xf=BS_LPfilter(x,fc,fs)
% function xf=BS_LPfilter(x,fc,fs)
%
% Standard Low-pass filter
% Requires Signal Processing Toolbox

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

[b,a]= butter(6,fc*2/fs);    %  Low pass filter coefficients
xm= mean(x);
xf= filtfilt(b,a,x-xm)+xm;    %  Low-pass filtered trace

return