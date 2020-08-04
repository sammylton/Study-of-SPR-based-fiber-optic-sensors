% This code uses matrix manipulations to plot Electric field intensity vs 
% distance from prism for the Kretschmann configuration of SPR sensors. 
% There are multiple layers with different thickness. 
% Each layer has a different magnetic permeability and a different 
% permittivity and thus a different refractive index.
% This code contains the block to find the resonance angle. 
% I think that one can also use the metal thickness optimisation code block
% i.e. ensure min(Reflectivity) < 0.01.
% This is for normalised electric field and not normalized magnetic field 
% as said on Pg8, Shalabney paper. This can be changed by changing the 
% field calculation step to the commented one right below it in the code.
clear; % removes all variables from the workspace
tic % "tic...toc" tells runtime
global np lambda u % so that their values can be used in the functions too 
%__________________________________________________________________________
%                               INPUT DATA
%__________________________________________________________________________    
%DataEmailedByAkhileshSir02Feb2020
u=pi*4e-7; % Magnetic Permeability of non-magnetic substances (assumed for every layer)
lambda=633e-9; % Wavelength of incident light in vacuum
% Wavelength affects refractive index and optimal thickness to be chosen 
% Optimal thickness usually keeps min(Reflectivity)< 0.01
N=4; % Number of layers including analyte and excluding prism
np=1.732; % Refractive Index of Prism (j=0)
na=1.33; % Refractive Index of Analyte (j=N layer)
% Refractive Index array of all the layers after prism : 
n=[sqrt(-26.767+24.086i);sqrt(-17.7372+0.6723i);sqrt(14.7097+0.1879i);na];
d=1e-9.*[10;40;10]; % Thickness chosen of layers between prism and analyte
% since analyte is assumed to continue to infinity
resotheta=72.546; % found Resonance Angle % Remove if block to find added
% TO CHECK, PUT 
%{
Shalabney Fig4a (only top curve)
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=2;
d=43e-9; but optimal distance is 46.98e-9(46nm) from FindMinThickness.m
n=[0.1325+4.0203*1i;na];%refractive indices of layers including that of analyte at wavelength=633nm
resotheta=54.6231; % 54.576 for d=46.98nm
%}
%{
% Shalabney Fig6b
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=3;
d=1e-9.*[43,10.5]; 
n=[0.1325+4.0203*1i;3.8354+0.0245*1i;na];
% resotheta=79.007; % resonance angle
resotheta =79.007; % 81 giving correct amplitude (as in Shalabney paper)
%}
% BLOCK FOR FINDING RESONANCE ANGLE ANGLE CAN BE ADDED RIGHT HERE
z=0:0.5e-9:2*(sum(d)); % row for thickness and RI arrays are column 
%__________________________________________________________________________
% Block for generating separate Zj for each layer starts.
%__________________________________________________________________________
% separate z are created because we want field at different values of z.
% visualise, different 3D 2X2matrices for each layer side-by-side, 3rd dimension being z. 
% when multiplied for each z, they give the final matrix. 
% the matrix for layer 2 doesnt impact the multiplication result for the z's it is identity matrix at.
% At every z, the multiplication depicts the propagation matrix calculated upto that very z
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
%__________________________________________________________________________
% Block for generating separate Zj for each layer ends
%__________________________________________________________________________
%Initiating matrices so that they don't change dimension during the loop
%__________________________________________________________________________
q=zeros(1,N);
b=zeros(N,numel(z));
M=zeros(2,2,numel(z),N);%3D-z values, 4D-layers
M0=zeros(2,2,N);
%__________________________________________________________________________
for j=1:N % loop on layers
    q(j)=qj(resotheta,n(j));
    b(j,:)=bj(resotheta,n(j),Z(j,:));
    Mtemp=[cosd(b(j,:));((sind(b(j,:))).*(1i.*q(j)));((sind(b(j,:))).*(1i./q(j)));cosd(b(j,:))];
% Since reshape works column-wise starting from  the first matrix-layer, we need to swap row 2 and 3 of Mm1 matrix in order to get correct 2X2 propagation matrix
    Mtemp=reshape(Mtemp,2,2,numel(z));
    M(:,:,:,j)=Mtemp;
    M0(:,:,j)=M(:,:,numel(z),j);% same matrix-layer for all z, hence one less Dimension from M
    % Last 3D component of the 4D matrix gives the matrix to calculate reflectivity 
    % the last layer of this matrix is superfluous
end

% CALCUATION OF THE PRODUCT OF ALL PROPAGATION MATRICES
M0tot(:,:)=eye(2,2);
for j=1:N-1 % because analyte layer is not included in reflectivity calculation
    % fourth matrix-layer of M0 is not required
    M0temp=reshape(M0(:,:,j),2,2);
    M0tot(:,:)=(M0tot(:,:))*(inv(M0temp));%new matrix multiplies to the right
end
Mtemp1=zeros(2,2,numel(z));
Mtot(:,:,:)=M(:,:,:,1);
for k=1:numel(z)
    for j=2:N
%         Mtemp=reshape(M(:,:,k,j),2,2,numel(z));
        Mtemp1(:,:,k)=M(:,:,k,j);
        Mtot(:,:,k)=(Mtemp1(:,:,k))*(Mtot(:,:,k));%downwards the layer, its matrix multiplies to the left side
    end
end
qp=qj(resotheta,np); % Value of function qj for prism
r=((M0tot(1,1)+(M0tot(1,2)).*(q(N))).*qp-(M0tot(2,1)+(M0tot(2,2)).*(q(N))))./((M0tot(1,1)+(M0tot(1,2)).*(q(N))).*qp+(M0tot(2,1)+(M0tot(2,2)).*(q(N))));
Exjz=-((Mtot(2,1,:)).*((1+r)./qp)+(Mtot(2,2,:)).*(1-r));% field calculation
% Exjz=-((Mtot(2,1,:)).*(1+r)+(Mtot(2,2,:)).*((1-r).*qp));
Exz(:)=Exjz;% Exjz=reshape(Exjz,[1,numel(z)]);%both serve same
modExjz2=(abs(Exz)).^2;
figure(1)% to distinguish different figures just in case
plot(z,modExjz2,'.-');
%__________________________________________________________________________
%                       PLOT VLINES AT EVERY BOUNDARY 
%__________________________________________________________________________
hold on;
sumd=0; % Initializing
for j=1:N-1
    sumd=sumd+d(j);
    plot([sumd sumd],[0 max(modExjz2)],'m--')
    hold on;
end
hold off;
%__________________________________________________________________________
xlabel('Distance from prism interface (m)'); 
ylabel('|E_x|^2');
legend(strcat(num2str(resotheta),'(deg)'),'Location','northeast');
toc % "tic...toc" tells runtime

function bofj = bj(theta,nj,zj)
% BOFJ cosoftheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2))); % cosine of the angle light ray makes at each interface.
    global lambda np 
    bofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((360*nj*zj)./lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
end

function qofj = qj(theta,nj)
% QOFJ cosoftheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%cos of angle in metal
    global u np
    qofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((sqrt(u))/nj);%qj in analyte %used cos of angle in that layer
end 