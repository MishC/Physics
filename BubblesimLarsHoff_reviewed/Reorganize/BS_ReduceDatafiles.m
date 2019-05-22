function BS_ReduceDatafiles
% function BS_ReduceDatafiles
%
% Sort results of Bubblesim calculation
% Remove unnecessary fields, to save space file
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

addpath(cd);
oldpath=cd;

%--- Select directory ---
filetype= '*.mat';
[filename,pathname]= uigetfile(filetype,'Select a file in the directory' );
cd(pathname);
set (findobj('Tag', 'resultpath'), 'String', pathname );

files = dir(filetype);
Nfiles= length(files);
drawnow;

n=0; m=0;
UpdateWaitbar (0);
for k=1:Nfiles
  WriteMessage(sprintf('Removing fields: File %d of %d  %s', k, Nfiles, files(k).name ) )
  UpdateWaitbar (k/Nfiles);
  drawnow
  newres= load( files(k).name );

  %--- Remove unnecesary fields  ---
  if isfield(newres,'pulse')
    n=n+1;
    particle  = newres.particle;
    pulse     = newres.pulse;
    simulation= newres.simulation;

    pulse     = rmexistfield(pulse, {'t','p'}  );
    simulation= rmexistfield(simulation, {'p'} );

    for l=1:length(simulation)
      simulation(l).a= simulation(l).a(:,1);
    end

    save ( files(k).name, 'particle', 'pulse', 'simulation' );
  else
    m=m+1;
    WriteMessage(sprintf('Rejected: %s', files(k).name ) )
  end
end

cd(oldpath);
WriteMessage(sprintf('Finished. %d files total, %d reduced', Nfiles, n ) )

return


%--- WRITEMESSAGE: Overrides m-file function with same name ---
% $$$ function WriteMessage(message)
% $$$ disp(message)
% $$$ return
