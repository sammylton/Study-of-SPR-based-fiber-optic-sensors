clear;
tic
nlambda=1;
lambda=linspace(0,600e-9,nlambda)';
%nm=0.53+10.43i;%real part is normal refractive index and complex part
nm=0.12+3.75i;%real part is normal refractive index and complex part
%is the decay term
%nm=0.1325;
np=1.732;
dna=1e-5;
na1=1.325;
na2=na1+dna;
na=[na1;na2];
u=pi*4e-7;%magnetic permebility of non-magnetic substances(assumed)
dm=47.3e-9;%thickness of metal layer
ntheta=1000000;%number of values of theta
theta=linspace(44,60,ntheta);%44 to 55 for chosen lambdas
costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));
bdm=(costheta).*((360*nm*dm)./lambda);
qp=(cosd(theta)).*((sqrt(u))/np);
qm=costheta.*(sqrt(u)/nm);
%M=[cosd(bdm);-((sind(bdm)).*(1i./qm));-((sind(bdm)).*(1i.*qm));cosd(bdm)];
M0(:,:,1)=cosd(bdm);
M0(:,:,2)=-((sind(bdm)).*(1i./qm));
M0(:,:,3)=-((sind(bdm)).*(1i.*qm));
M0(:,:,4)=cosd(bdm);
qa=zeros(1,ntheta,2);
indexmin=zeros(nlambda,2);
r=zeros(nlambda,ntheta,2);
R=zeros(nlambda,ntheta,2);
resoangle=zeros(1,2);
for k=1:2
    qa(:,:,k)=(sqrt(1-((sind(theta)).^2).*((np^2)./((na(k)).^2)))).*((sqrt(u))./(na(k)));
    r(:,:,k)=((M0(:,:,1)+(M0(:,:,2)).*(qa(:,:,k))).*qp-(M0(:,:,3)+(M0(:,:,4)).*(qa(:,:,k))))./((M0(:,:,1)+(M0(:,:,2)).*(qa(:,:,k))).*qp+(M0(:,:,3)+(M0(:,:,4)).*(qa(:,:,k))));
    R(:,:,k)=(abs(r(:,:,k))).^2;
    %plot(theta,R);
    %xlabel('Incidence angle (deg)');
    %ylabel('Reflectivity');
    %legend({'na=1.325','na=1.335'},'Location','southwest');
    for l=1:nlambda
        %Below is the block to find extrema
        indexmin(l,k) = find(min(R(l,:,k)) == R(l,:,k));
        resoangle(l,k) = theta(indexmin(l,k));
        %Above is the block to find extrema
    end
end
S=abs((resoangle(:,2)-resoangle(:,1))./dna);
figure(1)
plot(lambda,S,'.');
%plot(theta,R(:,:,1));% check
toc



