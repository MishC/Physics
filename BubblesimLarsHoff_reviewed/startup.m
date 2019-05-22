% Startup-file for 'Bubblesim' package
%
% 1) Define globals
% 2) Add program directories to Matlab's search path

% Lars Hoff, NTNU, Dept. of Telecommunications
% Trondheim, Norway

global BubblesimPath

BubblesimPath= cd;   % Global containing root path for program package

%==== DEFINE PATHS =========================================
if not(exist('BubblesimPath'))
  BubblesimPath= 'C:\Bubblesim';
end

% Add program directories to to Matlab's search path
addpath ( BubblesimPath,                            '-begin' );
addpath ( sprintf('%s\\Utilities', BubblesimPath ), '-begin' );
addpath ( sprintf('%s\\ODE',       BubblesimPath ), '-begin' );

% Start Bubblesim GUI menu 
BS_BubbleMenu