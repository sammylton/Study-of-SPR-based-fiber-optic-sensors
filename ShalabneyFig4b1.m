%For normalised magnetic field. 
clear;
tic
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
nm=0.1325+4.0203*1i;
np=1.732;
na=1.33;
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
dm=43e-9;%thickness of metal layer
k=0;
%{
%block to remove 1e9 starts
while 1 
    k=k+1;
    if eq(rem(dm*(10^k),1),0)
         break       
    end
end
%block to remove 1e9 ends
% but the above block does not work for 50e-9
%}
z=0:1e-9:5*dm;

%made it certain that dm lies in z vector
%z vector has to contain dm
%metal is from z=0 to dm
za=z-dm; za(za<=0)=0;%< or <= make no difference
%M2 is the reason of the discontinuity and the fact that it appears after z=dm(M2 acts like Identity before and at z=dm)
zdm=z; zdm(z>0)=dm; %zdm(1:1,1:m)=dm; is not so general as these two lines
zm=z; zm(z>dm)=dm;%> or >= make no difference
theta=[54.6231;55.6231;57.6231;53.6231;51.6231];
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
M=((M2(:,:,3)).*(M1(:,:,1))+(M2(:,:,4)).*(M1(:,:,3))).*((1+r)./qp)+((M2(:,:,3)).*(M1(:,:,2))+(M2(:,:,4)).*(M1(:,:,4))).*(1-r);%For normalised magnetic field.
%Ex=-M;

Mm=M;%Ex=-Mm
%Mm(:,(n+1:m))=0;
Mm(:,z>dm)=0;
Ma=M;
%Ma(:,(1:n))=0;%less general as n can be changed up in the code
Ma(:,z<=dm)=0;
%Either Mm or Ma has to be non zero at z=dm, otherwise you see no points at z=dm the figure
% Conceptually Ma should be zero at z=dm
Ezm=Mm.*(np/nm).*(sind(theta)./costheta1);
Eza=Ma.*(np/na).*(sind(theta)./costheta2);
Ez=Ezm+Eza;
%modEx2=(abs(M)).^2;
modEz2=(abs(Ez)).^2;
znm=z.*(1e9);%for x-axis to be nm scale
dmnm=dm*1e9;%for x-axis to be nm scale 
plot(znm,modEz2(1,:),'b.-',znm,modEz2(2,:),'r.--',znm,modEz2(3,:),'g.--',znm,modEz2(4,:),'r.-',znm,modEz2(5,:),'g.-',[dmnm dmnm],[0 max(modEz2(1,:))],'m--');%find a general method to create a vertical line
% . are being plotted to scrutinise discontinuity 
xlabel('Distance from prism interface (nm)'); 
ylabel('|Ez|2');
legend({strcat(num2str(theta(1)),'(deg)'),strcat(num2str(theta(2)),'(deg)'),strcat(num2str(theta(3)),'(deg)'),strcat(num2str(theta(4)),'(deg)'),strcat(num2str(theta(5)),'(deg)')},'Location','northeast');
toc