%For normalised magnetic field, hence the wrong graph
clear;
tic
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
nm=0.1325+4.0203*1i;
np=1.732;
na=1.33;
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
dm=43e-9;%thickness of metal layer
n=dm*1e9;
m=5*n+1;%5 so as to match with the given graph
z=linspace(0,(m-1)*1e-9,m);%not writing simply 5*dm to make it general and certain that dm lies in z vector
%metal is from z=0 to dm
za=z-dm;
za(za<=0)=0;% M2 is like Identity before z=dm
zdm(1:1,1:m)=dm;
zm=z;
zm(z>=dm)=dm;%> or >= makes no difference
theta=[54.6233;55.6233;57.6233;53.6233;51.6233];
%theta=53.6233;
costheta1=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));%angle in metal
costheta2=sqrt(1-((sind(theta)).^2).*((np^2)/(na^2)));%angle in analyte
bdm=(costheta1).*((360*nm/lambda).*zdm);%to be used in calculation of r
bm=(costheta1).*((360*nm/lambda).*zm);%360->2*pi when cosd->cos, sind->sin
%bm(z>=dm)=bdm;% M1 becomes constant after z=dm 
ba=(costheta2).*((360*na/lambda).*za);
qp=(cosd(theta)).*((sqrt(u))/np);
qa=(costheta2).*((sqrt(u))/na);
qm=costheta1.*(sqrt(u)/nm);
%M0=[cosd(bdm);-((sind(bdm)).*(1i/qm));-((sind(bdm)).*(1i*qm));cosd(bdm)];% Matrix for calculation of r
M0(:,:,1)=cosd(bdm);
M0(:,:,2)=-((sind(bdm)).*(1i./qm));
M0(:,:,3)=-((sind(bdm)).*(1i.*qm));
M0(:,:,4)=cosd(bdm);
%M1=[cosd(bm);((sind(bm)).*(1i/qm));((sind(bm)).*(1i*qm));cosd(bm)];% Matrix for propagation in metal 
M1(:,:,1)=cosd(bm);
M1(:,:,2)=((sind(bm)).*(1i./qm));
M1(:,:,3)=((sind(bm)).*(1i.*qm));
M1(:,:,4)=cosd(bm);
%M2=[cosd(ba);((sind(ba)).*(1i/qa));((sind(ba)).*(1i*qa));cosd(ba)];% Matrix for propagation in analyte
M2(:,:,1)=cosd(ba);
M2(:,:,2)=((sind(ba)).*(1i./qa));
M2(:,:,3)=((sind(ba)).*(1i.*qa));
M2(:,:,4)=cosd(ba);
%r=((M0(1,:)+(M0(2,:)).*qa).*qp-(M0(3,:)+(M0(4,:)).*qa))./((M0(1,:)+(M0(2,:)).*qa).*qp+(M0(3,:)+(M0(4,:)).*qa));
%M=((M2(3,:)).*(M1(1,:))+(M2(4,:)).*(M1(3,:))).*(1+r)+((M2(3,:)).*(M1(2,:))+(M2(4,:)).*(M1(4,:))).*((1-r).*qp);
r=((M0(:,:,1)+(M0(:,:,2)).*qa).*qp-(M0(:,:,3)+(M0(:,:,4)).*qa))./((M0(:,:,1)+(M0(:,:,2)).*qa).*qp+(M0(:,:,3)+(M0(:,:,4)).*qa));
M=((M2(:,:,3)).*(M1(:,:,1))+(M2(:,:,4)).*(M1(:,:,3))).*(1+r)+((M2(:,:,3)).*(M1(:,:,2))+(M2(:,:,4)).*(M1(:,:,4))).*((1-r).*qp);%For normalised magnetic field.
%Ex=-M;
modEx2=(abs(M)).^2;
znm=z.*(1e9);%for x-axis to be nm scale
dmnm=dm*1e9;%for x-axis to be nm scale 
plot(znm,modEx2(1,:),'b',znm,modEx2(2,:),'r',znm,modEx2(3,:),'g',znm,modEx2(4,:),'r--',znm,modEx2(5,:),'g--',[dmnm dmnm],[0 3e-6],'m--');%find a general method to create a vertical line
%,[dmnm dmnm],[0 max(modEx2(1,:))
xlabel('Distance from prism interface (nm)'); 
ylabel('|Ex|2');
legend({'54.6233(deg)','55.6233(deg)','57.6233(deg)','53.6233(deg)','51.6233(deg)'},'Location','northeast');
toc