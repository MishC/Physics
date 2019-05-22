% Diagonalize 1D Schrödinger equation -1/2 (d/dz)^2 - 1/((z-p)^2+b) in a basis of box states (sin(n*pi*x/L) for large
% L. The 200 lowest states are written out including energies and couplings. Number of x points are chosen
%such that fastest oscillatory function has 10 points per oscillation
% and Number of basis functions are determined by a maximum highest energy emax 
clear all; close all;
emax=1.0;  
L=2e+3; N=sqrt(2*emax)*L/pi; N=floor(N) 
p=L/2; b=1.4; 
dx=2*L/(N*pi)/5; Nx=L/(dx)
x=linspace(0,L,Nx); n=1:N; xp=1./sqrt((x-p).^2+b);
Ham=zeros(N,N) + diag((n.*pi).^2/(2*L^2));
N
for j=1:N
if mod(j,30) == 0
 j
end
  for k=j:N
    Ham(j,k)=Ham(j,k) - 2/L*sum(sin(j*pi*x/L).*sin(k*pi*x/L).*xp)*dx;  Ham(k,j)=Ham(j,k); 
    Zij=2/L*sum(sin(j*pi*x/L).*sin(k*pi*x/L).*x)*dx;
    if abs(Zij) > 1e-10
      l=l+1;
      indx(l)=j; indy(l)=k; Z(l)=Zij;
    end; 
  end;
end;
[V,D]=eig(Ham); E=diag(D); 
min(E)
% p-gauge: 
%for j=1:l
%  Z(j)=(E(indx(j))-E(indy(j)))*Z(j);
%end
save E3.dat E -ascii
ut=[indx' indy' Z'];
n=size(ut); n=n(1)
fid = fopen('Z3.dat','w'); 
for j=1:n
  fprintf(fid,'%5i  %5i  %14.9e\n',ut(j,:));
end
fclose(fid)
%save Z.dat ut -ascii
% Create the 50 lowest wavefunctions:
psi=zeros(50,Nx);
for j=1:50
  for k=1:N
    psi(j,:) = psi(j,:) + V(k,j).*sqrt(2/L).*sin(k*pi*x/L);
  end
[j, E(j) sum(abs(psi(j,:)).^2)*dx]
end
% Interpolate them onto a suitable FFT mesh.
Nxfft=2^12; psifft=zeros(100,Nxfft);
xfft=linspace(0,L,Nxfft); dxfft=xfft(2)-xfft(1);
for j=1:50
  psifft(j,:)=spline(x,psi(j,:),xfft);
  [j, E(j) sum(abs(psifft(j,:)).^2)*dxfft]
end
% save wavefunction at FFT points:
format long
cd in
save xfft3.dat xfft -ascii -double
save wfft3.dat psifft -ascii -double
cd ..

 

