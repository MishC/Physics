% Compute <x(t)>:

cd inF915
nbound=30; 
load wf.dat; psib=wf(1:nbound,:); load xf.dat; x=xf; dx=x(2)-x(1);
x=x-mean(x); % to have things centered in the middle
cd ../outF915

j=5;  s=char(j+48); stR = strcat('e1',s,'_Vhigh_psiR.dat'); stI = strcat('e1',s,'_Vhigh_psiI.dat');
psiR=load(stR); psiI=load(stI); clear i; psi=psiR+i*psiI;
load e15_High_T.dat; T=e15_High_T; T=T(:,1); T=T'
nT=size(psi); nT=nT(1); n

% plot expectation value of full wavefunction:
X_med=zeros(1,nT);
for k=1:nT
  X_med(k)= sum(abs(psi(k,:)).^2.*x)*dx;
  
end
figure(1); hold on; plot(T,X_med,'r')

% plot expectation value of continuum part of the wavefunction:

X_medC=zeros(1,nT);
for k=1:nT
  Psi_B=0;
  for j=1:nbound  
    Psi_B = Psi_B+sum(psib(j,:).*psi(k,:)*dx)*psib(j,:);
  end
  Psi_C=psi(k,:)-Psi_B;

  X_medC(k)= sum(abs(Psi_C).^2.*x)*dx;
end
figure(1); hold on; plot(T,X_medC,'b')

%figure
%for k=1:26
%  subplot(10,3,k); plot(x,abs(psi(k,:)).^2)
%end

