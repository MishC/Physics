function s= BS_WriteIndex(filename)
% function s= BS_WriteIndex(filename)
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% NO-7491, Trondheim, Norway
% January 2001


Mega = 1e6;
micro= 1e-6;

if filename ==0
  filetype= '*.mat';
  [filename,pathname]= uigetfile(filetype,'Select index file' );
  cd(pathname);
end

load (filename);

N= length(index);

for k= 1:N
  s{k,1}= index(k).filename;
  s{k,2}= index(k).particle.a0/micro;
  s{k,3}= index(k).pulse.A/Mega;
  s{k,4}= index(k).pulse.Nc;
  s{k,5}= index(k).pulse.f0/Mega;
  s{k,6}= index(k).pulse.envelope.name;
end

return