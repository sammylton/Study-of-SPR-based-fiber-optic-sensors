%For normalised electric field.
clear;
tic
global np lambda na u
% Values chosen, in order of incresing probability of getting changed

u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=4;%number of layers including analyte and excluding prism(as given in Shalabney paper)
d=1e-9.*[10;40;10]; 
n=[-26.767+24.086i;-17.7372+0.6723i;14.7097+0.1879i;na];%refractive indices of layers including that of analyte at wavelength=633nm
theta=87.2130; %found to be the angle of max reflectivity form the program ShalabneyFig3a4.m

% TO CHECK, PUT 
%{
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=2;
d=43e-9; 
n=[0.1325+4.0203*1i;na];%refractive indices of layers including that of analyte at wavelength=633nm
theta=54.6231;
%}
%{
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.332;
N=3;
d=1e-9.*[160,55]; 
n=[1.5968+0.02231i;0.13461+3.9882i;na];%refractive indices of layers including that of analyte at wavelength=633nm
theta=68.876;
%}
z=0:1e-9:2*(sum(d));

% Block for generating separate Zj for each layer starts.
Z=zeros(N,numel(z));
for j=1:N % Z of layer 1 starts from the 2nd row index.
    sumd=0;
    for k=1:j-1
        sumd=sumd+d(k);
    end
    Z(j,:)=z-sumd;
    ztemp=Z(j,:);%ztemp has been introduced for the next line (condition) can't do with double set of parentheses.
    ztemp(ztemp<0)=0;
    if j==N 
        Z(j,:)=ztemp;
        break % for the analyte has no limits
    end
    ztemp(ztemp>d(j))=d(j);
    Z(j,:)=ztemp;
end
% Block for generating separate Zj for each layer ends

qp=qj(theta,np);

q=zeros(N,numel(z));
b=zeros(N,numel(z));
M=zeros(2,2,numel(z),N);%can be initialised from all elements being 0
M0=zeros(2,2,N);

for j=1:N
    q(j,:)=qj(theta,n(j));
    b(j,:)=bj(theta,n(j),Z(j,:));
    Mtemp=[cosd(b(j,:));((sind(b(j,:))).*(1i.*q(j,:)));((sind(b(j,:))).*(1i./q(j,:)));cosd(b(j,:))];
% Since reshape works column-wise starting from  the first matrix-layer, we need to swap row 2 and 3 of Mm1 matrix in order to get correct 2X2 propagation matrix
    Mtemp=reshape(Mtemp,2,2,numel(z));% Important step
    M(:,:,:,j)=Mtemp;%Mtemp was created because M(:,:,:,j) is changing sizes which reshape didn't accept
    % M(:,:,:,j)'s is the propagation matrix of respective layer j.
    M0(:,:,j)=M(:,:,numel(z),j);% same matrix-layer for all z, hence one less Dimension from MMtemp
    % M0(:,:,j) is the inverse of complete z propagation matrix of layer j(z = end of layer j)
end
%{ 
M0tot(:,:)=eye(2,2);
for j=1:N-1
    % fourth matrix-layer of M0 is not required
    M0temp=reshape(inv(M0(:,:,j)),2,2);
    M0tot(:,:)=M0tot(:,:)*M0temp; % new multiplication should occur to the right
end
%}

Mtot=M(:,:,:,1);% for the nested loop starts from j=2
%{
Mtemp=zeros(2,2);
for k=1:numel(z)
    for j=2:N
        Mtemp=reshape(M(:,:,k,j),2,2);
%         Mtemp(:,:,k)=M(:,:,k,j);
        Mtot(:,:,k)=Mtemp*reshape(Mtot(:,:,k,1),2,2);%downwards the layer, its matrix multiplies to the left side
    end
end
%}

% below is a check loop for originally emailed problem, also for check values aforementioned in the code
M0tot=inv(M0(:,:,1))*inv(M0(:,:,2))*inv(M0(:,:,3));
for k=1:numel(z)
    Mtot(:,:,k)=reshape(M(:,:,k,4),2,2)*reshape(M(:,:,k,3),2,2)*reshape(M(:,:,k,2),2,2)*reshape(M(:,:,k,1),2,2);
end
%{
M0tot=inv(M0(:,:,1));
for k=1:numel(z)
    Mtot(:,:,k)=reshape(M(:,:,k,2),2,2)*reshape(M(:,:,k,1),2,2);
end
%}
%{
M0tot=inv(M0(:,:,1))*inv(M0(:,:,2));
for k=1:numel(z)
    Mtot(:,:,k)=reshape(M(:,:,k,3),2,2)*reshape(M(:,:,k,2),2,2)*reshape(M(:,:,k,1),2,2);
end
%}

r=((M0tot(1,1)+(M0tot(1,2)).*q(4)).*qp-(M0tot(2,1)+(M0tot(2,2)).*q(4)))./((M0tot(1,1)+(M0tot(1,2)).*q(4)).*qp+(M0tot(2,1)+(M0tot(2,2)).*q(4)));
%M=((M2(3,:)).*(M1(1,:))+(M2(4,:)).*(M1(3,:))).*(1+r)+((M2(3,:)).*(M1(2,:))+(M2(4,:)).*(M1(4,:))).*((1-r).*qp);
Exjz=-((Mtot(2,1,:)).*((1+r)./qp)+(Mtot(2,2,:)).*(1-r));
Exz(:)=Exjz;% Exjz=reshape(Exjz,[1,numel(z)]);%both serve same
modExjz2=(abs(Exz)).^2;
znm=z.*(1e9);%for x-axis to be nm scale
% dbynm=d.*1e9;%for x-axis to be nm scale 
figure(1)
plot(znm,modExjz2,'.-');
hold on;
count=0;
for j=1:N-1
    sumd=sumd+d(j)*1e9;
    plot([sumd sumd],[0 max(modExjz2)],'m--')
    hold on;
    count=count+1;
end
hold off;
%plot(znm,modEx2(1,:),'b',znm,modEx2(2,:),'r',znm,modEx2(3,:),'g',znm,modEx2(4,:),'r--',znm,modEx2(5,:),'g--',[dmnm dmnm],[0 25],'m--');%find a general method to create a vertical line
%,[dmnm dmnm],[0 max(modEx2(1,:))],'m--'
xlabel('Distance from prism interface (nm)'); 
ylabel('|Ex|2');
% legend(strcat('(nm1)=',num2str(n(1))),strcat(num2str(n(2)),'(nm2)'),strcat(num2str(n(3)),'(nd)'),strcat(num2str(n(4)),'(na)'),strcat(num2str(lambda),'(nm)'),strcat(num2str(theta),'(deg)'),'Location','northeast');
legend(strcat(num2str(theta),'(deg)'),'Location','northeast');
% legend(strcat(num2str(lambda),'(nm)'),'Location','northeast');
toc
%{
function cos=costheta(ri,thetaj)
    global np
    cos=sqrt(1-((sind(thetaj))^2)*((np^2)/(ri^2)));
end
%}

function bofj = bj(theta,nj,zj)
% BJ cosoftheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%cos of angle in metal
    global lambda np 
    bofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((360*nj*zj)./lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
end

function qofj = qj(theta,nj)
% QJ cosoftheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%cos of angle in metal
    global u np
    qofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((sqrt(u))/nj);%qj in analyte %used cos of angle in that layer
end 