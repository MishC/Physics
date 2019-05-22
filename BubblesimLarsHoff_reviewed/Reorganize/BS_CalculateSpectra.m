function BS_CalculateSpectra
% function BS_CalculateSpectra
%
% Calculate spectra of all files in directory and save to separate file
%

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

addpath(cd);
oldpath=cd;
SpectrumDir= 'Spectra';

%--- Select directory ---
filetype= '*.mat';
[filename,pathname]= uigetfile(filetype,'Select a file in the directory' );
addpath(pathname)

%--- Find mat-files ---
cd(pathname)
files = dir(filetype);
Nfiles= length(files);
drawnow;

set (findobj('Tag', 'resultpath'), 'String', cd );

%--- Loop thrpough all files ---
n=0; m=0;
UpdateWaitbar (0);

for k=1:Nfiles
  WriteMessage(sprintf('Fourier transforming: File %d of %d  %s', k, Nfiles, files(k).name ) )
  UpdateWaitbar (k/Nfiles);
  drawnow
  newres= load( files(k).name );

  if isfield(newres,'simulation')
    particle  = newres.particle;
    pulse     = newres.pulse;
    simulation= newres.simulation;

    n=n+1;
    for l=1:length(simulation)
      if isfield(simulation(l),'pr')
        Nfft= 2^nextpow2( length(simulation(l).pr) ); 
	Fpr = fft(simulation(l).pr, Nfft );

        simulation(l).Fpr = Fpr;
        simulation(l).Nfft= Nfft;

	%--- Plot result for inspection ---
	testplot= 1;
	if testplot
	  figure(3);
	  f = [0:Nfft-1]/Nfft*simulation(l).fs;
          Pr= 20*log10(abs(Fpr));  Pr= Pr-max(Pr);
       
	  plot(f, Pr, 'b.-')
          axis([0 20e6 -60 0])
	  xlabel('Frequency [Hz]')
	  ylabel(['Spectral power [dB]'])
	  hold off
	  drawnow
	  %keyboard
	end
      else
	WriteMessage(sprintf('Contents not recognized: %s', files(k).name ) )
      end
    end
    save ( files(k).name, 'particle', 'pulse', 'simulation' );
  end
end

cd(oldpath);
WriteMessage(sprintf('Finished. %d files total, %d reduced', Nfiles, n ) )

return


%--- WRITEMESSAGE: Overrides m-file function with same name ---
% $$$ function WriteMessage(message)
% $$$ disp(message)
% $$$ return
