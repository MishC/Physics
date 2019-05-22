function BS_UpdateWaitbar (x)
% function BS_UpdateWaitbar (x)
%
% Draw and update wiat-bar. x= 0...1


f = findobj('Tag', 'waitbar');
p = findobj(f,'Type','patch');
l = findobj(f,'Type','line');

if x==0
  subplot(f);
  cla

  xpatch = [ 0 0 0 0 ];
  ypatch = [ 0 0 1 1 ];
  xline  = [ 1 0 0 1 1];
  yline  = [ 0 0 1 1 0];

  p = patch(xpatch,ypatch,'r','EdgeColor','r','EraseMode','normal' ); % Clear old 
  set (p, 'EraseMode','none' )
  l = line(xline,yline,'EraseMode','none', 'color', 'k');

else

  if isempty(f) | isempty(p) | isempty(l), 
    error('Couldn''t find waitbar handles.'); 
  end

  xpatch = get(p,'XData');
  xpatch = [xpatch(2) x x xpatch(2)];
  set(p,'XData',xpatch')

  xline = get(l,'XData');
  set(l,'XData',xline);
end


return
