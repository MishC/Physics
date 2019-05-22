% Compute 10 ionization probabilities e0-e9 as P_ion = 1-P_norm-P_exc
% As there is a 1-2% norm loss in these very long calculations min(P_ion) is set manually to 0
clear all; close all;
cd in
nbound=20; nT=21;
load wfft8.dat; psib=wfft8(1:nbound,:); load xfft8.dat; x=xfft8; dx=x(2)-x(1);
Probs=zeros(10,3);
cd ../out15_A;

% open files and compute probs in a loop:
for j=0:9
   
  s=char(j+48); stR = strcat('e',s,'_psiR.dat'); stI = strcat('e',s,'_psiI.dat'); 
  psiR=load(stR); psiI=load(stI); clear i; psi=psiR(nT,:)+i*psiI(nT,:);
  norm=sum(abs(psi).^2)*dx; 
  if norm > 1 
    norm=1;
  end;
  Probs(j+1,1)=1-norm;
  
% start calculate bound state prob. for this particular field strength:
  Pe=0; 
    
  for k=1:nbound  
    Pe=Pe+abs(sum(psib(k,:).*psi)*dx).^2;
  end
  if Pe > 1 
    Pe=1;
  end;
  Probs(j+1,2)=Pe; Probs(j+1,3)=(1-Pe);

end;
Probs
cd .. 




