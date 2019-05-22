% Computes yearly production of Energy from Pelamis
% Power data are loaded from power_pelamis.dat.
% Wave and wind data are loaded from either wavedat_new.dat 
% FGN 11.09.06
% Last revision 06.11.14 

clear all;
close all;
%clc;

% INPUT
% Specifies which part of teh wave data file is to be used
Start_line= 1;                             % First data line to use in analyses (remove heading of input file)  
End_line = -1;                             % Last data line to use in analyses (<0: All data)   
T_1T02 = 1.20;                             % Ratio between T-01 and T02 for Jonswap (approximate values)
T_1Tp  = 0.85;                              % Ratio between T-01 and Tp for Jonswap  (approximate values)
% Read power characteristics
power_pelamis = load('power_pelamis.dat');  % Contains power curves for Pelamis. 
                                            % Note these data are for a 
                                            % 1250 kW unit. Performs a simple
                                            % scaling to 750 kW, see below
[Mp,Np] = size(power_pelamis);
% Read wave data
winddat = load('wavedat_new.dat');          % Wind and wave data file 

% End input
%***********************************************

% Organize the power data
Inp_Hs = power_pelamis(2:Mp,1);             % Siginficant wave height
Inp_T01 = power_pelamis(1,2:Np);            % Energy wave period T-01
Inp_Tp = Inp_T01/T_1Tp;                     % Converting the power table from energy period to Tp values (approximate)
Inp_power = power_pelamis(2:Mp,2:Np);       % Reading power values
Inp_power = Inp_power*750/1250;             % Scaling of capacity to 750kW

% Plot power data
figure
surf(Inp_Hs,Inp_Tp,Inp_power')
shading('interp')
title('Input Power Characteristic versus sign wave height and peak period')
ylabel('T_{p} (sec)')
xlabel('H_s (m)')
zlabel('P (kW)')
hold on

if End_line<0 
    End_line = size(winddat,1);
end
Nline0= End_line-Start_line+1;

% Assign wave data data 
month = winddat(Start_line:End_line,1);
day   = winddat(Start_line:End_line,2);
clock = winddat(Start_line:End_line,3);
minute = winddat(Start_line:End_line,4);
Hm0 = winddat(Start_line:End_line,5);
Tp = winddat(Start_line:End_line,6)';
% U = winddat(Start_line:End_line,8);
% Udir= winddat(Start_line:End_line,9);

%remove invalid numbers
KK = find(Hm0>=0 & Tp'>=0);    % Keeps samples where the wave height is non-negative.
month = month(KK);
clock = clock(KK);
day = day(KK);
Hm0n = Hm0(KK); 
Tpn = Tp(KK);

End_line = length(Hm0n);
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
tekst = ['Number of lines skipped (Hm0<0): ', num2str(DNline)];
disp(tekst)

Start_year = year(1);       
End_year = year(Nline);         

% Computes production for all measured wave conditions. Use 2D linear 
% interpolation in power table :
Power= zeros(Nline,1);          % Assign matrix 
Power_wav = Power;
for k=1:Nline
Power(k) = interp2(Inp_Tp,Inp_Hs,Inp_power,Tpn(k),Hm0n(k),'linear',0);        % Values outside the given interval are set to zero
end
Power_wav= ((1/(64*pi)*1025*9.80665^2*1.20).*Tpn'.*Hm0n.^2); %Power in waves

% Plots all production data on the power curve.
% plot3(Hm0,Tz,Power,'k+')


% Production per month.
mo=1;
while mo<= 12
    KK = find(month==mo);
    pow1 = Power(KK);      % power for given month,ignores NaN
    pow_wav= Power_wav(KK); 
    Hm0m = Hm0n(KK);             % Wind velocity for given month, ignores NaN   
    N_smon(mo) = size(pow1,1);
    Pow1_mon(mo)= mean(pow1);
    Pow_wav_mon(mo)=mean(pow_wav);
    Hm0_mon(mo) = mean(Hm0m);
    Hm0_max_mon(mo) = max(Hm0m);
    Hm0_min_mon(mo) = min(Hm0m);
    Hm0_std_mon(mo) = std(Hm0m);
    mo=mo+1;
end

%Capacity factors
Rated =  max(max(Inp_power))                % maximum power production
Cap1_mon = Pow1_mon./Rated                  % Capacity factor based upon averaging monthly data
Cap1_tot = mean(Power)./Rated               % Capacity factor based upon averaging all data
                    
        
MONTH = linspace(1,12,12);
figure
plot(MONTH,Cap1_mon,'Linewidth',2);
xlabel('Month');
ylabel('Capacity factor');
title ('Capacity factor per month');
legend('Based upon T_p')

figure
plot(MONTH,Pow_wav_mon/1000,'Linewidth',2);
xlabel('Month');
ylabel('Energy flux (kW/m)');
title ('Energy flux in waves per month');

% plot([1 12] ,[Cap_tot Cap_tot],'b--')
% hold on
figure
plot(MONTH,Hm0_mon,'Linewidth',2);
xlabel('Month');
ylabel('Mean significant wave height (m)');
title ('Mean Wave height per month');

figure
[haxes,hline1,hline2] =plotyy(MONTH,Hm0_max_mon,MONTH,Hm0_min_mon,'plot','plot');
set(hline1,'linewidth',2)
set(hline2,'linewidth',2)
xlabel('Month');
axes(haxes(1));
ylabel('H_{m0 max} (m)');
axes(haxes(2));
ylabel('H_{m0 min} (m)');
title ('Max and min significant wave height per month');

figure
plot(Hm0n,Tpn,'b+')
xlabel('Hmo, (m)')
ylabel('Tp, (sec)')
title ('Tp versus Hm0 for all wave conditions used');
axis([ 0 inf 0 inf])

figure
plot(Hm0n,Power_wav,'r+')
xlabel('Hmo, (m)')
ylabel('J, (kW/m)')
title ('Power versus Hm0 for all wave conditions used');
axis([ 0 inf 0 inf])


