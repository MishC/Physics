function filename = BS_WriteMessage(message, varargin)
% function filename = BS_WriteMessage(message, varargin)
%
% Simulate oscillation of bubble in acoustic field
%
% Write message in 'message' field
% 
% If 'messagetype' is present, the message is written to 
% the specified field
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

if isempty(varargin)
  field = 'message';
else
  field=varargin{1};
end

h= findobj('Tag', field);
if not(isempty(h))
  set(h, 'String', message )
end

return