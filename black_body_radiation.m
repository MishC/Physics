c0=2.997*10.^8;% m/s speed of light in vaccum
h=6.625*10.^-34;% J.s Planck constant 
k=1.38*10.^-23;%  T/K Boltzmann constant
n=1;%refravtive index of the medium .
step=  2.8429e-007;
lambda=[10.^-7:step:2*10.^-5];
T=1;%  Temperature initiatlisation
for i=300:50:800
    Mi=(2*pi.*h.*(c0.^2))./((n.^2).*(lambda.^5).*(exp((h.*c0)./(n.*k.*T.*lambda.*i))-1));
plot(lambda,Mi)
hold on
grid
end
xlabel('lambda')
ylabel('M(lambda)')
title('Monochromatic emittance')
legend('T=300','T=350','T=400','T=450','T=500','T=550','T=600','T=650','T=700','T=750','T=800')
hold off