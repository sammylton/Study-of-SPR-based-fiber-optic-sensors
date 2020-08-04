% Most general code to find resonance angle
clear;
tic % "tic...toc" tells runtime
global lambda np u  

% Shalabney Fig6b data
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=3;
d=1e-9.*[43;10.5]; 
n=[0.1325+4.0203*1i;3.8354+0.0245*1i;na];%refractive indices of layers including that of analyte at wavelength=633nm

%__________________________________________________________________________
% BLOCK TO BE ADDED TO OTHER CODES STARTS HERE
%__________________________________________________________________________
theta=0:1e-3:90;%angle in prism(angle of incidence)
Mtot=zeros(2,2,numel(theta)); % initializing for the very next loop
for l=1:numel(theta) % initialization loop for total characteristic matrix
    %Preallocationg in this loop for the next loop
    Mtot(:,:,l)=eye(2,2);
end
q=zeros(N-1,numel(theta));
b=zeros(N-1,numel(theta));
for j=1:N-1
    q(j,:)=qj(theta,n(j));
    b(j,:)=bj(theta,n(j),d(j));
    Mjtemp=[cosd(b(j,:));-((sind(b(j,:))).*(1i.*q(j,:)));-((sind(b(j,:))).*(1i./q(j,:)));cosd(b(j,:))];
    Mj= reshape(Mjtemp,2,2,numel(theta));
    for m=1:numel(theta)
        Mtot(:,:,m)=(Mtot(:,:,m))*(Mj(:,:,m));
    end
end
qp=qj(theta,np);qpr=(reshape(qp,1,1,numel(theta)));   
qa=qj(theta,na);qar=(reshape(qa,1,1,numel(theta)));
% reshaped to suit .* and ./ operations with TCmatrix elements(1X1Xnumel(theta)) correctly in the calculation of r
r=(((Mtot(1,1,:))+(Mtot(1,2,:)).*qar).*qpr-((Mtot(2,1,:))+(Mtot(2,2,:)).*qar))./(((Mtot(1,1,:))+(Mtot(1,2,:)).*qar).*qpr+((Mtot(2,1,:))+(Mtot(2,2,:)).*qar));%reflection coefficient
r=reshape(r,1,numel(theta)); % 1X1Xnumel(theta) to 1Xnumel(theta) for plottable form
R=(abs(r)).^2;%reflectivity
%Below is the block to find extrema
indexmin = find(min(R) == R); % to find resonance angle index
RA = theta(indexmin); % to find resonance angle from index
%__________________________________________________________________________
% BLOCK TO BE ADDED TO OTHER CODES ENDS HERE
%__________________________________________________________________________
toc % "tic...toc" tells runtime
%__________________________________________________________________________
% FUNCTION DECLARATION
%__________________________________________________________________________
% costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%nj varies in layers
function bofj = bj(theta,nj,dj)
    global lambda np 
    bofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((360*nj*dj)./lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
end
function qofj = qj(theta,nj)
    global u np
    qofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((sqrt(u))/nj);%qj in analyte %used cos of angle in that layer
end 