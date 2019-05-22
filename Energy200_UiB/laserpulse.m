function [y] = laserpulse(t)
w=7.5991e-5;
A=3.8894e-5;
%t=0:400:42341;

y=-A.*sin(w.*t);

% if t>0
% y=-A.*sin(2*w.*t)
% end
end