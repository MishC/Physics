function xs=BS_Scale125(x)
% function xs=BS_Scale125(x)
%
% Set axis scale on a 1-2-5 sequence
%
%  x : Axis maximum
% xs : Vector of length 5 containing ticks
%                      

% Lars Hoff
% NTNU, Dept. of Telecomm.
% N-7491 Trondheim, Norway
% 7 Jan 2000


amax    = max(x);
amaxexp = floor( log10(amax) );
amaxmant= amax/(10^amaxexp);
ascale= [ 1.5  2.0  2.5  3.0  4  5  6  8  10 ];
aticks= [ 0.5  0.5  0.5  0.5  1  1  1  2   2 ];
i=find( amaxmant<ascale );

amax = ascale(min(i));
atick= aticks(min(i));

xs= [0:atick:amax]'*10^amaxexp;

return






