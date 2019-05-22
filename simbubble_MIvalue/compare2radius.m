clear all;clf,clc;
global f R0 P0 Pv Pa ki PI rou sigma mu omiga 

f=0.5*10^6;
mi=0.30;
PI=3.1415926;                      
rou=998;  
ki=1.0;
sigma=0.072;                        
mu=1.*10^(-3);                            
Pv=2340; 
P0=101350;
Pa=sqrt(f)*mi*1e3;
R0=0.5*10^(-6);
fs=1*10^7;
omiga=2*PI*f;
time=1./fs;
tf=100.*time;
%solve
options=odeset('RelTol',1e-8,'AbsTol',[1e-7 1e-7]);
[t,y]=ode15s('diffode01',[0 tf],[R0 0],options,'Stats','on');
p1=plot((t),y(:,1)/R0,'k');xlabel('t [s]');
ylabel('Relative radius R(t)/R0');
hold on;
R0=4*10^(-6);
options=odeset('RelTol',1e-7,'AbsTol',[1e-7 1e-7]);
[t,y]=ode15s('diffode01',[0 tf],[R0 0],options,'Stats','on');
p2=plot((t),y(:,1)/R0,'-r');legend('R0=3e-6')
MII=strcat('MI= ',num2str(mi));
legend('R0=0.5e-6','R0=3e-6');title(MII)