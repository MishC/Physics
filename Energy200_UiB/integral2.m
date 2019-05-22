w=7.5991e-5;
A=3.8894e-5;
 t=linspace(-42341,0,400);
t1=linspace(0,42341,400);
 yy=-A.*sin(w.*t);

yx=-A.*sin(2.*w.*t);
%end
intr=trapz(t,yy)
intr1=trapz(t1,yx)
format long
difference=intr/intr1

