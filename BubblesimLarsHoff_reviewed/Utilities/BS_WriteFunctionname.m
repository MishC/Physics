function filename = BS_WriteFunctionname
% function filename = BS_WriteFunctionname
%
% Simulate oscillation of bubble in acoustic field
%
% Write name of currently active function
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

[stack n]= dbstack;
caller= stack(2).name;

h= findobj('Tag', 'functionname');
if not(isempty(h))
  set(h, 'String', caller  )
end

return



