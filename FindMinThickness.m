% Find explanation in ElectricFieldIntensityVsDistanceFromPrism.m code.
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
dm=40e-9:1e-11:50e-9;dm=dm';%thickness (so that Rmin<0.01)of metal layer for wavelength=632nm(given in Shalabney paper)
theta=54:1e-4:55;%angle in prism(angle of incidence)
costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));%cos of angle in metal
bdm=(costheta).*(((360*nm).*dm)./lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
qp=(cosd(theta)).*((sqrt(u))/np);%qj in prism (qj used in Shalabney paper)
qa=(sqrt(1-((sind(theta)).^2).*((np^2)/(na^2)))).*((sqrt(u))/na);%qj in analyte %used cos of angle in analyte
qm=costheta.*(sqrt(u)/nm);%qj in metal
%M=[cosd(b);-((sind(b)).*(1i./qm));-((sind(b)).*(1i.*qm));cosd(b)];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
M(:,:,1)=cosd(bdm);
M(:,:,2)=-((sind(bdm)).*(1i./qm));%qm and qdm are same because qj does not contain a z in its expression
M(:,:,3)=-((sind(bdm)).*(1i.*qm));
M(:,:,4)=cosd(bdm);
%r=((M(1,:)+(M(2,:)).*qa).*qp-(M(3,:)+(M(4,:)).*qa))./((M(1,:)+(M(2,:)).*qa).*qp+(M(3,:)+(M(4,:)).*qa));%reflection coefficient 
r=((M(:,:,1)+(M(:,:,2)).*qa).*qp-(M(:,:,3)+(M(:,:,4)).*qa))./((M(:,:,1)+(M(:,:,2)).*qa).*qp+(M(:,:,3)+(M(:,:,4)).*qa));
R=(abs(r)).^2;%reflectivity
%{
plot(theta,R(1,:));
xlabel('Incidence angle (deg)') 
ylabel('Reflectivity') 
%}
%{
%Below is the block to find extrema
indexmin = find(min(R(k,:)) == R(k,:)); 
thetamin = theta(indexmin); 
Rmin = R(indexmin);
indexmax = find(max(R) == R);
thetamax = theta(indexmax); 
Rmax = R(indexmax); 
%}
%{
R=R';
[Rmin,indexmin]=min(R);
%Rmin=Rmin';
indexmin2=find(min(Rmin)==Rmin);
Rminmin=Rmin(indexmin2);
dmmin=dm(indexmin2);
thetamin=theta(indexmin(indexmin2));
plot(theta,R(:,indexmin2))
%}
indexmin = min(min(R));
[x,y]=find(R==indexmin);
Rmin=R(x,y);
dmmin=dm(x);
thetamin=theta(y);
plot(theta,R(x,:))

toc