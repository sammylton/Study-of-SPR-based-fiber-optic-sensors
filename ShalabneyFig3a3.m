clear;
tic %tic (at the start) and toc (at the end) finds the runtime 
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
%changing wavelength to 1500nm from 633nm (no other changes) changes the whole nature of the graph for the refractive index of the metal changes
% Also, a different thickness is required to get minimum reflectivity (Rmin) under 0.01
% Changing thickness dm changes the resonance angle
nm=0.1325+4.0203i;%refractive index of metal at wavelength=633nm
%nm=0.53+10.43i;%refractive index of metal at wavelength=1500nm
%nm=0.12+3.75i;%refractive index of metal at wavelength=600nm
np=1.732;%refractive index of prism 
na=1.33;%refractive index of analyte
u=pi*4e-7;%magnetic permeability of non-magnetic substances(assumed)
%dm=43e-9;%thickness of metal layer for wavelength=633nm
%dm=47e-9;%thickness (so that Rmin<0.01)of metal layer for wavelength=632nm(given in Shalabney paper)
dm=46.98e-9;%thickness (so that Rmin<0.01, acc. to me it is making Rmin nearly vanish)of metal layer for wavelength=633nm approx.
%dm=48e-9;%thickness (so that Rmin<0.01)of metal layer for wavelength=600nm (approximated through the help of paper)
% dm=43.25e-9;%thickness of metal layer for wavelength=850nm(given in Shalabney paper)
% dm=29.75e-9;%thickness of metal layer for wavelength=1550nm(given in Shalabney paper)
%dm=30.69e-9;%thickness (so that Rmin<0.01, acc. to me it is making Rmin nearly vanish)of metal layer for wavelength=1500nm (approximated myself, not 10nm, not 15nm)
theta=0:1e-1:90;%angle in prism(angle of incidence)
costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));%cos of angle in metal
b=(costheta).*((360*nm*dm)/lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
qp=(cosd(theta)).*((sqrt(u))/np);%qj in prism (qj used in Shalabney paper)
qa=(sqrt(1-((sind(theta)).^2).*((np^2)/(na^2)))).*((sqrt(u))/na);%qj in analyte %used cos of angle in analyte
qm=costheta.*(sqrt(u)/nm);%qj in metal
M=[cosd(b);-((sind(b)).*(1i./qm));-((sind(b)).*(1i.*qm));cosd(b)];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
r=((M(1,:)+(M(2,:)).*qa).*qp-(M(3,:)+(M(4,:)).*qa))./((M(1,:)+(M(2,:)).*qa).*qp+(M(3,:)+(M(4,:)).*qa));%reflection coefficient 
R=(abs(r)).^2;%reflectivity
plot(theta,R);
xlabel('Incidence angle (deg)') 
ylabel('Reflectivity') 

%Below is the block to find extrema
indexmin = find(min(R) == R); 
thetamin = theta(indexmin); 
Rmin = R(indexmin);
indexmax = find(max(R) == R);
thetamax = theta(indexmax); 
Rmax = R(indexmax); 

toc