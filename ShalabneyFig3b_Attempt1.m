%This program deals with two values of na 
clear;
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
nm=0.1325+4.0203i;%real part is normal refractive index and complex part
%is the decay term
%nm=0.1325;
np=1.732;%refractive index of prism 
na=[1.325;1.335];%Make column vector of refractive index of analyte so that matrix dimension get 1Xmany(theta takes many values) to 2Xmany
u=pi*4e-7;%magnetic permeability of non-magnetic substances(assumed)
dm=43e-9;%thickness of metal layer
theta=linspace(0,90,1000);%angle in prism(angle of incidence)
costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));
b=(costheta).*((360*nm*dm)/lambda);
qp=(cosd(theta)).*((sqrt(u))/np);
qa=(sqrt(1-((sind(theta)).^2).*((np^2)./(na.^2)))).*((sqrt(u))./na);
qm=costheta.*(sqrt(u)/nm);
M=[cosd(b);-((sind(b)).*(1i./qm));-((sind(b)).*(1i.*qm));cosd(b)];
r=((M(1,:)+(M(2,:)).*qa).*qp-(M(3,:)+(M(4,:)).*qa))./((M(1,:)+(M(2,:)).*qa).*qp+(M(3,:)+(M(4,:)).*qa));
%r=((M11+M12.*qn).*q0-(M21+M22.*qn))./((M11+M12.*qn).*q0+(M21+M22.*qn));
R=(abs(r)).^2;
plot(theta,R);
xlabel('Incidence angle (deg)'); 
ylabel('Reflectivity');
legend({'na=1.325','na=1.335'},'Location','southwest');
%{
%Below is the block to find extrema
indexmin = find(min(R) == R); 
thetamin = theta(indexmin); 
Rmin = R(indexmin);
indexmax = find(max(R) == R);
thetamax = theta(indexmax); 
Rmax = R(indexmax); 
%Above is the block to find extrema
%}