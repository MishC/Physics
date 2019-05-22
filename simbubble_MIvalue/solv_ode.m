function p1=solv_ode(d)
global  f R0 P0 Pv Pa ki rou sigma mu omiga time step
                     
rou=998;  
ki=1.4;%ratio of heat capacities
sigma=0.073;                        
mu=1.*10^(-3); %Pa s                           
Pv=2340; 
P0=1.0*10^5;
MI=d;
%Pa=4*10^4;
R0=2.5*10^(-6);
f=.5*10^6;
omiga=2*pi*f;
Pa=sqrt(f/1e6)*MI*1e6;
time=1./f;%how much samples per sec
step=20;
ttime=step.*time;%number of steps
options=odeset('RelTol',1e-5,'AbsTol',[1e-5 1e-5], 'Stats', 'on');
[t,y]=ode15s('diffode01',[0 ttime],[R0 0.],options);

p1=plot((t*f),y(:,1)/R0);
xlabel('Time/T');ylabel('Solution R(t)/R0');
legend('R=R(t)/R0');
% figure(2)
% plot((t*f),p);
% 
% xlabel('Time/T');ylabel('Driving pressure');
