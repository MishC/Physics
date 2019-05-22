function y= BS_rmexistfield( var, field)
% function y= BS_rmexistfield( var, field)
%
% function y= rmexistfield( var, field)
%
% Remove field of struct-variable if it exists

% Lars Hoff, NTNU, Dept. of Telecom, NO-7491 Trondheim, Norway
% 02 Jan 2001

y= var;
for k=1:length(field)
  if isfield( y, field{k} )
    y= rmfield( y, field{k} );
  end
end

return
  