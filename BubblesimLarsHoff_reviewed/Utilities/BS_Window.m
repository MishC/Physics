function w = BS_Window( envelope, N, r )
% function w = BS_Window( envelope, N, r )
%
% Window of length N points
%
% envelope : string
%            'Rectangular'
%            'Triangular'
%	     'Hanning'
%            'Hamming'
%	     'Blackman'
%            'Costapered rectangular'
%
%        N : No. of points
%        r : Additional parameter used by some windows
%

%  Lars Hoff, NTNU, Dept. of Telecommunications 
%  N-7491 Trondheim, Norway
%  Oct. 2000
%  Modified Jan 2001 LH

n= [1:N]';              % Indexes         1,2,3 ... N
t= (n-1)./(N-1)-1/2;    % Relative scale  -0.5 ... +0.5
x= 2*pi*t;              % Radians         -pi ... pi

switch lower(envelope(1:6))
 case 'rectan', w = ones(N,1);
 case 'triang'; w = 1-2*abs(t);
 case 'hannin'; w = 0.50+0.50*cos(x); 
 case 'hammin'; w = 0.54+0.46*cos(x);
 case 'blackm'; w = 0.42+0.50*cos(x)+0.08*cos(2*x);

 case 'costap'; if (r>1)
                  error('Parameter r must be less than 1');
		  return
                end
		c = 1/2*(1-r);
		k2= find(t> c);
		t2= 1/2*(t(k2)-c)./(1/2-c);
		x2= 2*pi*t2;
                k1= N+1-k2;

		w = ones(N,1);
		w(k2)= 0.50+0.50*cos(x2);
		w(k1)= w(k2);

 otherwise, error ('Unknown window name')
end

return

