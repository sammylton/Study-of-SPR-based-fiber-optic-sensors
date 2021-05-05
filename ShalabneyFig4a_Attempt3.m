%For normalised electric field. In ShalabneyFig4a2, I used already calculated resonance angle
% And this code contains resonance angle finding block
% One can also use the metal thickness optimisation code block here
clear;
tic
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
%same wavelength used in the 0nm legend of Fig6b(read theory to know)
%nd=3.8354+0.0245i at lambda=633e-9
nm=0.1325+4.0203*1i;
np=1.732;
na=1.33;
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances

% One can also use the metal thickness optimisation code block here
dm=43e-9;%thickness of metal layer
%dm=46.98e-9;%thickness of metal layer

z=0:1e-10:5*dm;%for x-axis of required graph
%metal is from z=0 to dm

% To make M2 like Identity before z=dm
za=z-dm;
za(za<=0)=0;

zm=z;
zm(z>=dm)=dm;%I don't know how it can be z>=dm and should not compulsarily be zm>=dm 
%, maybe MATLAB has it inbuilt to check it as a condiition
%> or >= makes no difference


%block to find resonance theta for given values starts
theta_inc=linspace(0,90,1000000);%angle in prism(angle of incidence)
costheta=sqrt(1-((sind(theta_inc)).^2).*((np^2)/(nm^2)));%cos of angle in metal
b=(costheta).*((360*nm*dm)/lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
qp=(cosd(theta_inc)).*((sqrt(u))/np);%qj in prism (qj used in Shalabney paper)
qa=(sqrt(1-((sind(theta_inc)).^2).*((np^2)/(na^2)))).*((sqrt(u))/na);%qj in analyte %used cos of angle in analyte
qm=costheta.*(sqrt(u)/nm);%qj in metal
M=[cosd(b);-((sind(b)).*(1i./qm));-((sind(b)).*(1i.*qm));cosd(b)];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
r=((M(1,:)+(M(2,:)).*qa).*qp-(M(3,:)+(M(4,:)).*qa))./((M(1,:)+(M(2,:)).*qa).*qp+(M(3,:)+(M(4,:)).*qa));%reflection coefficient 
R=(abs(r)).^2;%reflectivity
indexmin=find(min(R) == R);
thetamin = theta_inc(indexmin); 
%block to find resonance theta for given values ends

%theta=[54.6231;55.6231;57.6231;53.6231;51.6231];
%theta=53.6233;
%In ShalabneyFig4a2, I used already calculated resonance angle
theta=[thetamin;thetamin+1;thetamin+3;thetamin-1;thetamin-3];

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
M=((M2(:,:,3)).*(M1(:,:,1))+(M2(:,:,4)).*(M1(:,:,3))).*((1+r)./qp)+((M2(:,:,3)).*(M1(:,:,2))+(M2(:,:,4)).*(M1(:,:,4))).*(1-r);%For normalised electric field.
%Ex=-M;
modEx2=(abs(M)).^2;
znm=z.*(1e9);%for x-axis to be nm scale
dmnm=dm*1e9;%for x-axis to be nm scale 
plot(znm,modEx2(1,:),'b',znm,modEx2(2,:),'r',znm,modEx2(3,:),'g',znm,modEx2(4,:),'r--',znm,modEx2(5,:),'g--',[dmnm dmnm],[0 25],'m--');%find a general method to create a vertical line
%,[dmnm dmnm],[0 max(modEx2(1,:))],'m--'
xlabel('Distance from prism interface (nm)'); 
ylabel('|Ex|2');
legend({strcat(num2str(theta(1)),'(deg)'),strcat(num2str(theta(2)),'(deg)'),strcat(num2str(theta(3)),'(deg)'),strcat(num2str(theta(4)),'(deg)'),strcat(num2str(theta(5)),'(deg)')},'Location','northeast');

toc