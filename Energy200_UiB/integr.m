clear;close
w=2.5333e-5;
A=3.8894e-5;
 %t=linspace(-42341,42341,400);


int=quadl('laserpulse', -42341,0, 1e-6,1);
int1=quadl('laserpuls', 0, 42341,1e-6,1);
int1=abs(int1);
format long
differenc=int/int1;

x=linspace(-42341,0,400);
x1=linspace(0,20000,400);
p1=laserpulse(x);
p2=laserpuls(x1);
for a=1:length(x)
iy(a)=quad('laserpulse',-42341,0);

iyy(a)=quad('laserpuls',0,42341/2);
end
figure(1)
plot(x,p1,'-r',x1,p2,'-r')
%hold on
% plot(x,iy,'g',x1,iyy,'g')


