function BS_ResampleFiles
% function BS_ResampleFiles
%
% Reasample all files in directory to constant rate
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
  WriteMessage(sprintf('Resampling: File %d of %d  %s', k, Nfiles, files(k).name ) )
  UpdateWaitbar (k/Nfiles);
  drawnow
  newres= load( files(k).name );

  if isfield(newres,'simulation')
    particle  = newres.particle;
    pulse     = newres.pulse;
    simulation= newres.simulation;

    n=n+1;
    for m=1:length(simulation)
      if isfield(simulation(m),'p')
	fs = 20*pulse(m).f0;
	[tr,pr,fs]= ConstantSampleRate( simulation(m).t, simulation(m).p, fs );
	simulation(m).fs= fs;
	simulation(m).pr= pr;
      
	%--- Plot result for inspection ---
	testplot= 0;
	if testplot
	  figure(3);
	  t= [0:length(simulation(m).pr)-1]/fs;

	  plot(simulation(m).t, simulation(m).p, 'b.-')
	  hold on
	  plot(t, simulation(m).pr, 'r.-')
	  hold off
	  drawnow
	  %keyboard
	end
	save ( files(k).name, 'particle', 'pulse', 'simulation' );
      else
	WriteMessage(sprintf('Already resampled: %s', files(k).name ) )
      end
    end
    %simulation= rmfield(simulation, {'p'}  );
  end
end

cd(oldpath);
WriteMessage(sprintf('Finished. %d files total, %d reduced', Nfiles, n ) )

return


%--- WRITEMESSAGE: Overrides m-file function with same name ---
% $$$ function WriteMessage(message)
% $$$ disp(message)
% $$$ return
