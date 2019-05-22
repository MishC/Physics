% Diagonalize 1D Schrödinger equation -1/2 (d/dz)^2 - 1/((z-p)^2+b) in a basis of box states (sin(n*pi*x/L) for large
% L. The 200 lowest states are written out including energies and couplings. Number of x points are chosen
%such that fastest oscillatory function has 10 points per oscillation
% and Number of basis functions are determined by a maximum highest energy emax 
clear all; %close all;
emax=6.0;  
L=127.75; N=sqrt(2*emax)*L/pi; N=floor(N) 
p=L/2; b=2.0; 
dx=2*L/(N*pi)/10; Nx=floor(L/(dx))
x=linspace(0,L,Nx); n=1:N; xp=1./sqrt((x-p).^2+b);
Ham=zeros(N,N) + diag((n.*pi).^2/(2*L^2));
N

for j=1:N
if mod(j,30) == 0
 j
end
  for k=j:N
    Ham(j,k)=Ham(j,k) - 2/L*sum(sin(j*pi*x/L).*sin(k*pi*x/L).*xp)*dx;  Ham(k,j)=Ham(j,k); 
  end;
end;
[V,D]=eig(Ham); E=diag(D); 
min(E)

save E3.dat E -ascii
% Create the 10 lowest wavefunctions:
psi=zeros(10,Nx);
for j=1:10
  for k=1:N
    psi(j,:) = psi(j,:) + V(k,j).*sqrt(2/L).*sin(k*pi*x/L);
  end
[j, E(j) sum(abs(psi(j,:)).^2)*dx]
end
% Interpolate them onto a suitable FFT mesh.
Nxfft=2^9; psifft=zeros(10,Nxfft);
xfft=linspace(0,L,Nxfft); dxfft=xfft(2)-xfft(1);
for j=1:10
  psifft(j,:)=spline(x,psi(j,:),xfft);
  [j, E(j) sum(abs(psifft(j,:)).^2)*dxfft]
end
% save wavefunction at FFT points:
format long
save xfft4.dat xfft -ascii -double
save wfft4.dat psifft -ascii -double


 

