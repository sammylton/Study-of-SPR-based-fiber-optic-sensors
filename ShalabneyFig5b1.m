%This is using the final code for Fig4a. For normalised electric field. The vector theta has been enlarged to find the resonance angle from the block added.
clear;
tic
lambda=632e-9;%wavelength of incident light in vacuum at which nm is as follows
nm=0.1325+4.0203*1i;
np=1.732;
na=1.33;
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
dm=46.98e-9;%thickness of metal layer
z=0:1e-9:7*dm;
%metal is from z=0 to dm

% To make M2 like Identity before z=dm
za=z-dm;
za(za<=0)=0;

zm=z;
zm(z>=dm)=dm;%I don't know how it can be z>=dm and should not compulsarily be zm>=dm %> or >= makes no difference

theta=linspace(54,55,10001)';%I know that it is going to come between these values and upto 4 digit after the decimals' precision
%theta=53.6233;
costheta1=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));%angle in metal
costheta2=sqrt(1-((sind(theta)).^2).*((np^2)/(na^2)));%angle in analyte
bdm=(costheta1).*((360*nm/lambda).*dm);%to be used in calculation of r
bm=(costheta1).*((360*nm/lambda).*zm);%360->2*pi when cosd->cos, sind->sin
%bm(z>=dm)=bdm;% M1 becomes constant after z=dm 
ba=(costheta2).*((360*na/lambda).*za);
qp=(cosd(theta)).*((sqrt(u))/np);
qa=(costheta2).*((sqrt(u))/na);
qm=costheta1.*(sqrt(u)/nm);
%M0=[cosd(bdm);-((sind(bdm)).*(1i/qm));-((sind(bdm)).*(1i*qm));cosd(bdm)];
% Matrix for calculation of r
M0(:,:,1)=cosd(bdm);
M0(:,:,2)=-((sind(bdm)).*(1i./qm));%qm and qdm are same because qj does not contain a z in its expression
M0(:,:,3)=-((sind(bdm)).*(1i.*qm));
M0(:,:,4)=cosd(bdm);
%M1=[cosd(bm);((sind(bm)).*(1i/qm));((sind(bm)).*(1i*qm));cosd(bm)];
% Matrix for propagation in metal 
M1(:,:,1)=cosd(bm);
M1(:,:,2)=((sind(bm)).*(1i./qm));
M1(:,:,3)=((sind(bm)).*(1i.*qm));
M1(:,:,4)=cosd(bm);
%M2=[cosd(ba);((sind(ba)).*(1i/qa));((sind(ba)).*(1i*qa));cosd(ba)];
% Matrix for propagation in analyte
M2(:,:,1)=cosd(ba);
M2(:,:,2)=((sind(ba)).*(1i./qa));
M2(:,:,3)=((sind(ba)).*(1i.*qa));
M2(:,:,4)=cosd(ba);
%r=((M0(1,:)+(M0(2,:)).*qa).*qp-(M0(3,:)+(M0(4,:)).*qa))./((M0(1,:)+(M0(2,:)).*qa).*qp+(M0(3,:)+(M0(4,:)).*qa));
%M=((M2(3,:)).*(M1(1,:))+(M2(4,:)).*(M1(3,:))).*(1+r)+((M2(3,:)).*(M1(2,:))+(M2(4,:)).*(M1(4,:))).*((1-r).*qp);
r=((M0(:,:,1)+(M0(:,:,2)).*qa).*qp-(M0(:,:,3)+(M0(:,:,4)).*qa))./((M0(:,:,1)+(M0(:,:,2)).*qa).*qp+(M0(:,:,3)+(M0(:,:,4)).*qa));

%block to find resonance angle starts
R=(abs(r)).^2;%reflectivity
indexmin = find(min(R) == R); 
thetamin = theta(indexmin); 
Rmin = R(indexmin);
%block to find resonance angle ends

M=((M2(:,:,3)).*(M1(:,:,1))+(M2(:,:,4)).*(M1(:,:,3))).*((1+r)./qp)+((M2(:,:,3)).*(M1(:,:,2))+(M2(:,:,4)).*(M1(:,:,4))).*(1-r);%For normalised electric field.
%Ex=-M;
modEx2=(abs(M)).^2;
znm=z.*(1e9);%for x-axis to be nm scale
dmnm=dm*1e9;%for x-axis to be nm scale 
%plot(znm,modEx2(1,:),'b',znm,modEx2(2,:),'r',znm,modEx2(3,:),'g',znm,modEx2(4,:),'r--',znm,modEx2(5,:),'g--',[dmnm dmnm],[0 25],'m--');%find a general method to create a vertical line
%,[dmnm dmnm],[0 max(modEx2(1,:))],'m--'
plot(znm,modEx2(indexmin,:));
xlabel('Distance from prism interface (nm)'); 
ylabel('|Ex|2');
%legend({strcat(num2str(theta(1)),'(deg)'),strcat(num2str(theta(2)),'(deg)'),strcat(num2str(theta(3)),'(deg)'),strcat(num2str(theta(4)),'(deg)'),strcat(num2str(theta(5)),'(deg)')},'Location','northeast');
toc