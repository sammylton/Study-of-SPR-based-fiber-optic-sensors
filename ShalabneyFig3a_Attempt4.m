% Most general code  for Fig3a
clear;
tic %tic (at the start) and toc (at the end) finds the runtime 
global lambda np u  
%{
% Values chosen, in order of increasing probability of getting changed
%DataEmailedByAkhileshSir02Feb2020
u=pi*4e-7;%magnetic permeability of non-magnetic substances(assumed for every layer)
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;%refractive index of prism 
na=1.33;%refractive index of analyte
N=4;%number of layers including analyte and excluding prism(as given in Shalabney paper)
n=[sqrt(-26.767+24.086i);sqrt(-17.7372+0.6723i);sqrt(14.7097+0.1879i)];%refractive index of metal at wavelength=633nm
d=1e-9.*[10;40;10];%thickness (so that Rmin<0.01, acc. to me it is making Rmin nearly vanish)of metal layer for wavelength=633nm approx.
%}
% TO CHECK, PUT
% Shalabney fig4a data 
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;%refractive index of prism
na=1.33;%refractive index of analyte
N=2;%number of layers including analyte and excluding prism
n=0.1325+4.0203i;%refractive index of metal at wavelength=633nm
% d = 46.98e-9;%thickness (so that Rmin<0.01)of metal layer
d = 43e-9; % thickness to be taken according to the paper
%{
% Shalabney Fig6b
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=3;
d=1e-9.*[43,10.5]; 
n=[0.1325+4.0203*1i;3.8354+0.0245*1i;na];%refractive indices of layers including that of analyte at wavelength=633nm
%}
theta=0:1e-4:90;%angle in prism(angle of incidence)
M=zeros(2,2,numel(theta));

for l=1:numel(theta)
    %Preallocationg M in this loop for the next loop
    M(:,:,l)=eye(2,2);
end

for j=1:N-1
    Mk=[cosd(bj(theta,n(j),d(j)));-((sind(bj(theta,n(j),d(j)))).*(1i.*qj(theta,n(j))));-((sind(bj(theta,n(j),d(j)))).*(1i./qj(theta,n(j))));cosd(bj(theta,n(j),d(j)))];
    Mj= reshape(Mk,2,2,numel(theta));
    for m=1:numel(theta)
        M(:,:,m)=(M(:,:,m))*(Mj(:,:,m));
    end
end
M11(:)=M(1,1,:);
M12(:)=M(1,2,:);
M21(:)=M(2,1,:);
M22(:)=M(2,2,:);
r=((M11+(M12).*(qj(theta,na))).*(qj(theta,np))-(M21+(M22).*(qj(theta,na))))./((M11+(M12).*(qj(theta,na))).*(qj(theta,np))+(M21+(M22).*(qj(theta,na))));%reflection coefficient
R=(abs(r)).^2; % Reflectivity
plot(theta,R,'.');
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
% costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%nj varies in layers 
function bofj = bj(theta,nj,dj)
    global lambda np 
    bofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((360*nj*dj)./lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
end
function qofj = qj(theta,nj)
    global u np
    qofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((sqrt(u))/nj);%qj in analyte %used cos of angle in that layer
end 