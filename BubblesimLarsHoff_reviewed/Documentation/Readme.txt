'BUBBLESIM'. INSTALLING AND RUNNING
-------------------------------------
Please be aware that this is an experimantal program, intended for research. 
It is far from the the quality standards of a commercial product.
I attempt to keep this file up to date, but some changes may still not
have been included. 
Any comments are welcome!


FUNCTION OF THE PROGRAM
-----------------------
The program 'Bubblesim' simulates the response of a bubble exposed to
an ultrasound pulse. 
The main purpose of the program is to simulate contrast
agents for medical ultrasound imaging. It is aimed at micrometer-sized
bubbles exposed to Megahertz-frequency ultrasound.
Only the radial oscillation mode is included, i.e. the oscillation is
spherically symmetric.
Nonlinear terms are fully included. Nonlinear results are caculated 
using time-domain simulations.
A visco-elastic shell is included. Thickness and visco-elastic
parameters of the shell must be specified by the user. 


SYSTEM REQUIREMENTS
---------------------
The original program was written in Matlab ver. 5.2 and Windows NT 4.0.
It has over time been updated to new Matlab and Windows versions.
The current program was written and tested under Matlab ver 7 and Windows XP. 

Any computer running Matlab ver 5.2 or newer should be able to run the program.
The graphical input and output may need modifications on non-Windows
computers. 

A Pentium or equivalent processor is recommended, to obtain reasonable
run-time of the simulations. The memory requirements depend on the
size of the traces to be simulated.

Functions from the Signal Processing Toolbox are used for filtering
and resampling. This is not critical for the calculations. If the 
Signal Processing Toolbox is not found, the calculations will be
performed without filtering, and a warning message is issued. 

The graphical displays were written for a 1280x1024 pixels screen,
running Windows 'Large Fonts'. The program has also been tested and
works acceptable at screen resolutions 1024x768 and 800x600 pixels, at
both 'Large Fonts' and 'Small Fonts'. 
However, several different screens exist, and the visual appearance of
the graphical windows may look strange at some screen/font settings. 


INSTALLATION
----------------
Unzip the files to selected folder. 
Do not extract to a path name containing spaces, (e.g. 
'c:\Program Files\'), as this may cause problems with Matlab.  
Recommended is to unzip to the root directory, e.g. 'c:\'.
A new folder 'Bubblesim' will be created under the selected folder.


RUNNING
-----------------
1) Start Matlab

2) Change to directory containing simulation files
   e.g. 'cd c:\Bubblesim'

3) Running from graphical user interface

   a) Run startup-program 'startup.m'.
      This sets up paths and starts the programs.
      Alternatively, enter Matlab command 'Bubblesim'
   
   b) Enter particle and pulse parameters

   c) Commands
      Display pulse: Display selected pulse
      Calculate    : Simulate bubble response

   d) Options
     'Inverted pulse'          : Include inverted pulse in simulation
     'Plot linear calculations': Include results of a linearized 
                                 model, for comparison

4) Running from batch file
   Program 'BubblesimBatch.m'
   - Written for running batch jobs
   - Allows several bubble diameters, pulse amplitudes etc.
   - Run manually by changing the Matlab-code in 
     'BubblesimBatch.m', then calling the program from Matlab.
   - Results are stored in .mat-files in an open format. These can 
     be loaded for later viewing and plotting results, and for 
     calculations and comparisons.



RESULTS
------------------------------------------------------------------

Results are stored in Matlab 'struct'-variables
          pulse: Driving ultrasound pulse
       particle: Parameters of the contrast agent particle, or bubble
         linear: Results of linear calculation
     simulation: Results of nonlinear simulation
                 simulation(1): Original pulse
                 simulation(2): Inverted pulse, if selected
          graph: Plotting and saving parameters


CONTENTS OF RESULT-VARIABLES
------------------------------------------------------------------

pulse
      pulse(1): Original pulse data
      pulse(2): Data of inverted pulse, if selected

  Fields:
      envelope: Pulse envelope name and corresponding Matlab command
             A: Pulse amplitude               [Pa]
            Nc: Pulse length                  [No. of cycles]
            f0: Pulse center frequency        [Hz]
            fs: Pulse sample rate             [1/s]
        invert: Include inverted pulse        [Boolean]
             t: Time vector                   [s]
             p: Driving pressure pulse vector [Pa]
        source: Name of pulse source: file name or 'Synthetic'
------------------------------------------------------------------

particle
            p0: Internal equilibrium pressure         [Pa]
           rho: Density of the surrounding liquid     [kg/m3]
            eL: Viscosity of the surrounding liquid   [Pas]
             c: Speed of sound in the liquid          [m/s]
            Kg: Thermal conductivity of the gas       [W/(mK)]
            rg: Density of gas                        [kg/m3]    
            Cp: Heat capacity of the gas at p=const   [J/(kgK)]  
         gamma: Adiabatic constant of the gas, Cp/Cv  [1]
            a0: Particle radius                       [m]
            ds: Shell thickness                       [m]
            Gs: Shell shear modulus                   [Pa]
            es: Shell viscosity                       [Ns/m2]
----------------------------------------------------------------------

linear 
             t: Time vector                  [s]
             x: Relative radial displacement [1]
        a(:,1): Particle radius              [m] 
        a(:,2): Particle velocity            [m/s] 
             p: Scattered pressure           [Pa]
             f: Frequency vector             [Hz]
           Hxp: Transfer function: Pressure in -> Radial displacement
           Hpp: Transfer function: Pressure in -> Pressure out
           hxp: Impulse response : Pressure in -> Radial displacement
           hpp: Impulse response : Pressure in -> Pressure out
----------------------------------------------------------------------

simulation
    simulation(1): Original pulse data
    simulation(2): Data of inverted pulse, if selected

  Fields:
           solver: ODE solver. Name and Matlab function
            model: Bubble model. Name and Matlab function
   thermaldamping: Thermal model for the gas
  displayprogress: Plot progress during simulation          [Boolean]
            etime: Elapsed CPU time for simulation             [s] 
    -- Unevenly sampled vectors, as returned from ODE solver --
                t: Time vector.        Unevenly sampled        [s]
           a(:,1): Particle radius.    Unevenly sampled        [m] 
           a(:,2): Particle velocity.  Unevenly sampled        [m/s] 
                p: Scattered pressure. Unevenly sampled        [Pa]
    -- Vectors resampled to constant rate fs --
               tr: Time vector.        Constant sample rate fs [s]
               pr: Scattered pressure. Constant sample rate fs [Pa]
               fs: Actual sample rate of simulation            [1/s]

----------------------------------------------------------------------

graph
       plotlinear: Include linear calculations   [Boolean]
          include: Graphs to include in plot
            title: Simulation title
       resultfile: File to store simulation results into
             tmax: Max. time to plot             [s]
             fmax: Max. frequency to plot        [Hz]
           figure: Figure no. where results are plotted
         incoming: Graph for driving pulse(s)          {figure, graph}
        scattered: Graph for scattered pulse           {figure, graph}
           radius: Graph for bubble radius             {figure, graph}
         velocity: Graph for wall velocity             {figure, graph}
          spectra: Graph for power spectra             {figure, graph}
         transfer: Graph for linear transfer functions {figure, graph} 
           symbol: Line styles used in plots
----------------------------------------------------------------------



FILE NAMING CONVENTION
----------------------------------------------------------------------

  Filename: CSyyyymmdd_nnn.mat

  Symbol  Description      Values 
  --------------------------------------------------------------------
   CS     File code        'Contrast Simulation'
   yyyy   Year             2000, 2001, 2002, ...
   m      Month            01, 02, 03, ... , 12
   dd     Day of month     01, 02, 03, ... , 31
   nnn    Measurement no.  001, 002, 003 ... , 999
   mat    Extension        Matlab file format
  --------------------------------------------------------------------

  Example: The 15th simulation result on June 25, 2001, is
stored in a file named 'CS20010625_015.mat'.



