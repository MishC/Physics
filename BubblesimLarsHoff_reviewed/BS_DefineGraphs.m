function [graph] = BS_DefineGraphs(graph)
% function [graph] = BS_DefineGraphs(graph)
%
% Bubblesim: Simulate bubble response
% Define graphs and plotting parameters
%
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

BS_WriteFunctionname;

%=== DEFINE RESULT WINDOWS ===========================================
graph.figure= 2;

if not( ishandle(graph.figure) )
  figure(graph.figure);
  screensize = get (0, 'screensize');
  figsize = screensize(3:4)/2;
  set (graph.figure, ...
        'Name',                  'Results',  ... 
	'NumberTitle',           'off', ...
	'units',                 'points', ...
	'position',               [figsize(1)-5 figsize(2)-30 figsize ], ...
	'DefaultAxesBox',        'On', ...      
	'DefaultAxesFontName',   'Arial', ...
	'DefaultAxesFontWeight', 'normal', ...
	'DefaultAxesFontSize',    8, ...
	'DefaultLineLineWidth',   1        );
else
  figure(graph.figure); 
  clf;
end


%--- Determine which graphs to plot where ---
switch sum(graph.include);
 case  1,     M=1;  N=1;
 case  2,     M=1;  N=2;
 case {3,4},  M=2;  N=2;
 case {5,6},  M=2;  N=3;
end
  
diagram(1).code= 'incoming';    diagram(1).title= 'Driving Pulse';
diagram(2).code= 'scattered';   diagram(2).title= 'Scattered Pulse';
diagram(3).code= 'radius';      diagram(3).title= 'Bubble Radius';
diagram(4).code= 'velocity';    diagram(4).title= 'Bubble Wall Velocity';
diagram(5).code= 'spectra';     diagram(5).title= 'Power Spectra';
diagram(6).code= 'transfer';    diagram(6).title= 'Transfer Functions';

n=0;
for k=1:length(graph.include)
  if graph.include(k)
    n=n+1;
    w= { graph.figure, subplot(M,N,n) };
    title(diagram(k).title);
  else
    w= 0;
  end
  eval( sprintf('graph.%s = w;', diagram(k).code ));
end


%=== DEFINE PLOT SYMBOLS =========================================
graph.symbol.incoming  = 'b-';
graph.symbol.linear    = 'k:';
graph.symbol.attenuated= 'b--';
graph.symbol.result    = {'r-';'b-'};

return














