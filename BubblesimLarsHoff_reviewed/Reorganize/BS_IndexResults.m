function BS_IndexResults
% function BS_IndexResults
%
% Sort results of Bubblesim calculation
% Place all parameters into one index file
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

addpath(cd);
oldpath=cd;

%--- Select directory ---
filetype= '*.mat';
[filename,pathname]= uigetfile(filetype,'Select a file in the directory to index' );
cd(pathname);
set (findobj('Tag', 'resultpath'), 'String', pathname );

files = dir(filetype);
Nfiles= length(files);
drawnow;

resultfile= sprintf('%sindex', files(1).name(1:5) );

n=0; m=0;
UpdateWaitbar(0);
for k=1:Nfiles
  pulse     = [];
  particle  = [];
  simulation= [];

  WriteMessage(sprintf('Indexing: File %d of %d  %s', k, Nfiles, files(k).name ) )
  UpdateWaitbar (k/Nfiles);
  drawnow
  newres= load( files(k).name );
  
  %--- Store fields to index  ---
  if isfield(newres,'pulse')
    n=n+1;

    index(n).particle  = newres.particle;
    index(n).pulse     = newres.pulse(1);
    index(n).simulation= newres.simulation(1);
    index(n).filename  = files(k).name;
    index(n).date      = files(k).date;

    index(n).pulse     = rmexistfield( index(n).pulse,      {'t','p'}     );
    index(n).simulation= rmexistfield( index(n).simulation, {'t','a','p','tr','pr', 'Fpr'} );
 
  else
    m=m+1;
    WriteMessage(sprintf('Rejected: %s', files(k).name ) )
  end
end

WriteMessage(sprintf('Saving results to %s', resultfile ) )
save (resultfile, 'index' );

cd(oldpath);
WriteMessage(sprintf('Finished. %d files total, %d indexed', Nfiles, n ) )

return


