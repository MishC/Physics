function filename = BS_MakeFilename
% function filename = BS_MakeFilename
%
% Simulate oscillation of bubble in acoustic field
%
% Generate filename for results from 'Bubblesim' 
% Automatically generated from  date and measurement number 
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway


type  = 'CS';   % Measurement code for bubble simulation
time  = clock;
year  = mod(time(1),100);
month = time(2);
day   = time(3);

number  = 0;
occupied= 1;
while occupied
  number  = number+1;
  filename= sprintf('%2s%02d%02d%02d_%03d.mat', type, year, month, day, number );
  occupied= exist (filename, 'file');
end  
 
return