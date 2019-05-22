function [dy]=diffode01(t,y)
global  R0 P0 Pv Pa ki rou sigma mu omiga 
dy=zeros(2,1);
% p=drive(t);
% if p<=2*PI/omiga*p
%     p=Pa*exp(0.2)*sin(omiga*t);
% elseif 2*PI/omiga*p<p<8*PI/omiga*p
%     
%     p=Pa*sin(omiga*t);
% elseif 11>p>=8*PI/omiga*p
%         p=Nc/time*Pa*exp(-0.2)
% else p=0;
%     Nc/step*(1/omiga)*Pa*sin(omiga*t))
% end
Nc=5;
dy(1)=y(2);
dy(2)=((P0+2.*sigma./R0-Pv).*(R0./y(1)).^(3*ki)+Pv-2.*(sigma./y(1))-4*mu.*y(2)./y(1)-...
P0-Pa.*sin(omiga*t))./(rou.*y(1))-(1.5.*(y(2).^2)./y(1));

%y(1)=R
%y(2)= dR
%dy(2) =acceleration