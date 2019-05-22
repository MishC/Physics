% Computes yearly production of wind power from a wind turbine.
% Turbine power curve is located in teh file power_curves.dat.
% Wave and wind data are loaded from either wavedat_new.dat.
% FGN 11.09.06
% Last revision 06.11.14

clear all;
%close all;
%clc;
format compact

% INPUT
H  = 85;                                 % Heigth of turbine axis (m)
alfa = 0.133;                            % Heigth exp. U= U_10*(z/10)^alfa
Start_line= 1;                           % First data line to use in analyses (remove heading of input file)
End_line = -1;                           % Last data line to use in analyses (<0: All data)

%Read power curve
load power_curve.dat;                   % Contains power curve for wind turbine
                                        % This file contains data for a
                                        % 2.5MW turbine, scale it to 5 MW

% Wind Data
winddat = load('wavedat_new.dat');

% End input
%***********************************************
%
Inp_wi = power_curve(:,1);
Inp_po = power_curve(:,2)*2;           % Factor 2 to scale form 2.5 to 5 MW

scale = (H/10)^alfa;                % Correction of wind velocities from 10m to H

%**********************************************************************


if End_line<0
    End_line = size(winddat,1);
end
Nline0= End_line-Start_line+1;

% Assign data
month = winddat(Start_line:End_line,1);
day   = winddat(Start_line:End_line,2);
clock = winddat(Start_line:End_line,3);
minute = winddat(Start_line:End_line,4);
U = winddat(Start_line:End_line,8);
Udir= winddat(Start_line:End_line,9);

% Remove invalid numbers
% Keeps samples where the wind velocity is non-negative.
U_pos_ind = find(U >= 0);
month = month(U_pos_ind);
day = day(U_pos_ind);
clock = clock(U_pos_ind);
Un = U(U_pos_ind)*scale;           % Velocity in nacelle height

End_line = length(Un);
Nline= End_line-Start_line+1;
DNline = Nline0-Nline;
year= 1*ones(size(month));       %assigns all months to year 1
fclose all;

% write data interval used
tekst = [' Starting date and time:'];
disp (tekst)
tekst = [' Year - month - day - hour'];
disp (tekst)
date_1 = [year(1) month(1) day(1) clock(1)]
tekst = [' End date and time:'];
disp (tekst)
tekst = [' Year - month - day - hour'];
disp (tekst)
date_2 = [year(Nline) month(Nline) day(Nline) clock(Nline)]
tekst = ['Number of lines skipped (U<0): ', num2str(DNline)];
disp(tekst)

Start_year = year(1);
End_year = year(Nline);

% Computes production for all measured wind velocities, uses cubic interpolation:
Power = interp1(Inp_wi,Inp_po,Un,'cubic',0);        % Values outside the given interval are set to zero

%Computes distribution of Power production
hour=linspace(1,8760,size(Un,1));
figure
plot(hour,sort(Un),'Linewidth',3)
grid on
axis([0 8760 0 inf]);
xlabel('Hours')
ylabel('U_n (m/sec)')
title({'10 min average wind velocity below given level';'Hub heigth'})

figure
plot(hour,sort(Power),'Linewidth',3)
grid on
axis([0 8760 0 inf]);
xlabel('Hours')
ylabel('P (kW)')
title(['Number of hours with power production below given level'])

% Production per month.
mo=1;
while mo<= 12
    KK = find(month==mo);
    pow = Power(KK);            % power for given month,ignores NaN
    uw = Un(KK);                % Wind velocity for given month, ignores NaN
    N_smon(mo) = size(pow,1);
    Pow_mon(mo)= mean(pow);
    Un_mon(mo) = mean(uw);
    Un_max_mon(mo) = max(uw);
    Un_min_mon(mo) = min(uw);
    Un_std_mon(mo) = std(uw);
    mo=mo+1;
end

%Capacity factors
%Cap_year = Pow_year./max(Inp_po)
Cap_mon = Pow_mon./max(Inp_po)
Cap_tot = mean(Power)./max(Inp_po)
Mean_cap_month = mean(Cap_mon)                  % Differences in these mean values reflect different length of months and the effect of missing data
Mean_wind_tot = mean(Un)
Mean_wind_mon = mean(Un_mon)

% plot results

MONTH = linspace(1,12,12);
figure
[haxes,hline1,hline2] =plotyy(MONTH,Cap_mon,MONTH,Un_mon,'plot','plot','Linewidth',2);
set(hline1,'linewidth',2)
set(hline2,'linewidth',2)
xlabel('Month');
axes(haxes(1));
ylabel('Capacity factor');
axes(haxes(2));
ylabel('U_{mean} (m/sec)');

figure
[haxes,hline1,hline2] =plotyy(MONTH,Un_max_mon,MONTH,Un_min_mon,'plot','plot');
set(hline1,'linewidth',2)
set(hline2,'linewidth',2)
xlabel('Month');
axes(haxes(1));
ylabel('U_{max} (m/sec)');
axes(haxes(2));
ylabel('U_{min} (m/sec)');


