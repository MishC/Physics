function dy=RP(t,y)
global R0 P0 Pv Pa ki rou sigma mu omiga 
dy=zeros(2,1);
dy(1)=y(2);
dy(2)=((P0-Pv+2*sigma/R0)*(R0./y(1))^(3.*ki)-2*sigma/y(1)-4*mu*y(2)/y(1)+Pv-...
(P0+Pa*sin(omiga*t)))/(rou*y(1))-1.5*(y(2)^2)/y(1);
return

