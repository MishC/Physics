function y= BS_GetUserdata(h)
% function value= BS_GetUserdata(h)
%
% Read user data from GUI pop-up menu with handle 'h'.
% Returns 'UserData(value)', where 'value' is the selected pop-up 'value'

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway


n   = get ( h, 'value'    );  
data= get ( h, 'UserData' );
   
y= data(n);

return
