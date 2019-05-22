function [t,p]= BS_MakePulse(A,Nc,f0,fs,envelope)
% function [t,p]= BS_MakePulse(A,Nc,f0,fs,envelope)
%
%         A : Amplitude
%        Nc : No. of cycles
%        f0 : Centre frequency [Hz]
%        fs : Sample frequency [Hz]
%  envelope : Name of pulse envelope, windowing funtion
%
%   Construct ultrasound pulse
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

%-- Time ---
T = Nc/f0;               % Pulse length
dt= 1/fs;                % Sample interval
tp= (0:dt:T )';          % Time vector, containing oscillations
t = (0:dt:4*T )';        % Time vector, total

%--- Pulse ---
W  = BS_Window( envelope{1}, length(tp), envelope{2});
po = sin(2*pi*f0*tp) ;    % Carrier wave
po = A*po.*W;             % Pulse with envelope

%--- Place pulse ---
p= zeros( size(t) );
n= [1:length(po) ]; % Index to put oscillations into
p(n)= po;

return









