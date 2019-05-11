% cleaning
clear all; close all; clc;

%free bubbles in the water
%data 
global f R0 P0 Pv Pa ki PI rou sigma mu omiga MI
f=0.5*10^6;
MI=0.5;
PI=3.1415926;                      
rou=998;  
ki=1.0;
sigma=0.072;                        
mu=1.*10^(-3);                            
Pv=2340; 
P0=1.013*10^5;
Pa=sqrt(f/1e6)*MI*1e6;
R0=2.5*10^(-6);
fs=1*10^7;
omiga=2*PI*f;
time=1./fs;
tf=100.*time;
%solve
options=odeset('RelTol',1e-6,'AbsTol',[1e-8 1e-8],);
[t,y]=ode15s('RP',[0 tf],[R0 0],options,'Stats','on');
%plot
plot((t),y(:,1)/R0,'Linewidth',2);
p=num2str(Pa/1e3);mi=num2str(MI); nn= strcat('MI=', mi);
 rr=num2str(R0*1e6);
 pp1=[' and ',p,' kPa'];rr1=[' and ',rr,' um']
 
 titul=strcat(nn, pp1,rr1);,
xlabel('Time [s]');ylabel('R(t)/R0');title(titul);
hold on;
%plot([0,R0],[0 0],'bk', 'Linewidth',1);
plot(t,ones(size(t)),'k','Linewidth',1.5);hold on


