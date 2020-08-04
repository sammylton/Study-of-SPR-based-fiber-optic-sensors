% This code plots square of x-Electric field's amplitude vs distance from prism for the Kretschmann configuration of SPR sensors. 
% Having the paper(DOI: 10.1016/j.sna.2010.02.005) will help in deciphering the code for some symbols are matching.
% There are multiple layers in the sensor with different: 
% (1) Thickness (depends on the wavelength of incident light)
% (2) Magnetic Permeability 
% (3) Permittivity = (Refractive Index)^2 (changes with wavelength)
% Minimum reflectivity of incident light occurs at the resonance angle, given the thickness of each layer.
% This code may contain the block to find the resonance angle.
% I THINK that one can also do metal thickness optimization to ensure min(Reflectivity) < 0.01 by varying the resonance angle slightly around its found value for one thickness(nested loop on angle and thickness).
% Since the field depends on the (1)total propagation matrix, (2)reflection coefficient at the resonance angle, to find it at every z:
% (1) total propagation matrix (TPmatrix) at every z is needed: 
% Hence TPmatrix is of 2X2Xnumel(z) dimension(3D).
% TPmatrix(2X2Xnumel(z)) = Product of N Propagation matrices(2X2Xnumel(z)).
% Elements of each shell of a Propagation matrix(Pmatrix) depend upon properties of a layer and z values.
% Consider a 3D matrix for layer 2 as its propagation matrix.
% Its 2X2 "shells(2D)" become constant after z>(layer1+layer2) thickness and are Identity before z<(layer1) thickness.
% Hence separate 'Z' matrix(N X numel(z)) is created with a dimension added to the z vector. 
% Rows correspond to layers and the different Z rows are used for calculation of different Pmatrices of dimension 2X2Xnumel(z) for each layer.
% (2) To find the value of reflection coefficient at the resonance angle, the calculation of total characteristic matrix has to be done. 
% TCmatrix(2X2Xnumel(z)) = Product of N-1(except analyte) Characteristic matrices(2X2Xnumel(z)).
% Since, for each layer, Cmatrix=inverse(last 2D shell of Pmatrix), and all Z rows become constant at the z analyte starts, 
% TCmatrix=inverse(last shell of product of all Pmatrices except that of the analyte).
% 
% By Sameer Baheti under Dr. Akhilesh Kumar Mishra, IITR
% sbaheti@ph.iitr.ac.in
% dated
clear; % removes all variables from the workspace
tic % "tic...toc" tells runtime
global np lambda u % so that their values can be used in the functions too 
%__________________________________________________________________________
%                INPUT DATA FOR KRETSCHMANN CONFIGURATION
%__________________________________________________________________________    
%Data emailed by Dr. Akhilesh Kumar Mishra on 02Feb2020
u=pi*4e-7; % Magnetic Permeability of non-magnetic substances (assumed for every layer).
lambda=633e-9; % Wavelength of incident light in vacuum.
N=4; % Number of layers including analyte and excluding prism.
np=1.732; % Refractive Index of Prism (j=0).
na=1.33; % Refractive Index of Analyte (j=N layer).
% Refractive Index array of all the layers after prism: 
n=[sqrt(-26.767+24.086i);sqrt(-17.7372+0.6723i);sqrt(14.7097+0.1879i);na];
d=1e-9.*[10;40;10]; % Thickness chosen of layers excluding prism and analyte.
% Since analyte is assumed to continue to infinity, calculations are such that d array is one element shorter than refractive index array.
% RA=72.546; % Resonance Angle % Comment if block to find is added
% TO CHECK, PUT 
%{
% Shalabney Fig4a (only top curve)
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=2;
d=43e-9; but optimal distance is 46.98e-9(46nm) from FindMinThickness.m
n=[0.1325+4.0203*1i;na];%refractive indices of layers including that of analyte at wavelength=633nm
RA=54.6231; % 54.576 for d=46.98nm
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
RA=79.007; % but 81 giving correct amplitude (as in Shalabney paper)
%}
% BLOCK TO FIND RESONANCE ANGLE CAN BE ADDED RIGHT HERE.
%__________________________________________________________________________
% BLOCK TO FIND RESONANCE ANGLE STARTS HERE
%__________________________________________________________________________
theta=0:1e-3:90; % angle in prism (angle of incidence onto layer system)
Mtot=zeros(2,2,numel(theta));

% loop to initialize total characteristic matrix
for l=1:numel(theta)
    %Preallocationg M in this loop for the next loop
    Mtot(:,:,l)=eye(2,2);
end

for j=1:N-1
    %Mk=[cosd(bj(theta,n(j),d(j)));-((sind(bj(theta,n(j),d(j)))).*(1i./qj(theta,n(j))));-((sind(bj(theta,n(j),d(j)))).*(1i.*qj(theta,n(j))));cosd(bj(theta,n(j),d(j)))];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
    Mk=[cosd(bj(theta,n(j),d(j)));-((sind(bj(theta,n(j),d(j)))).*(1i.*qj(theta,n(j))));-((sind(bj(theta,n(j),d(j)))).*(1i./qj(theta,n(j))));cosd(bj(theta,n(j),d(j)))];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
    %Since reshape works column-wise starting from  the first matrix-layer,
    %we need to swap row 2 and 3 of Mk matrix in order to get correct 2X2 propagation matrix    
    Mj= reshape(Mk,2,2,numel(theta));
    for m=1:numel(theta)
        Mtot(:,:,m)=(Mtot(:,:,m))*(Mj(:,:,m));
    end
end
M11(:)=Mtot(1,1,:);
M12(:)=Mtot(1,2,:);
M21(:)=Mtot(2,1,:);
M22(:)=Mtot(2,2,:);
%reflection coefficient:
r=((M11+(M12).*(qj(theta,na))).*(qj(theta,np))-(M21+(M22).*(qj(theta,na))))./((M11+(M12).*(qj(theta,na))).*(qj(theta,np))+(M21+(M22).*(qj(theta,na))));
R=(abs(r)).^2;%reflectivity
indexmin = find(min(R) == R); 
RA = theta(indexmin); % Resonance Angle
%__________________________________________________________________________
% BLOCK TO FIND RESONANCE ANGLE ENDS HERE
%__________________________________________________________________________
z=0:0.5e-9:2*(sum(d));
%__________________________________________________________________________
% Block for generating separate Zj for each layer starts
%__________________________________________________________________________
% separate Z are created because we want field at different values of z.
% visualise, different 3D 2X2matrices for each layer side-by-side, 3rd dimension being z. 
% when multiplied for each z, they give the final matrix. 
% the matrix for layer 2 doesn't impact the multiplication result for the z's it is identity matrix at.
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
%__________________________________________________________________________
% TOTAL PROPAGATION MATRIX (PRODUCT OF ALL PROPAGATION MATRICES) CALCUATION
% TOTAL CHARACTERISTIC MATRIX CALCUATION
%__________________________________________________________________________
% Initializing matrices so that they don't change dimension during the loop
q=zeros(1,N);
b=zeros(N,numel(z));
M=zeros(2,2,numel(z),N); 
% 3D propagation matrices clubbed together in the fourth dimension to assist multiplication of matrices later in loop
for j=1:N % loop on layers
    q(j)=qj(RA,n(j)); % qj at resonance angle
    b(j,:)=bj(RA,n(j),Z(j,:)); % bj at resonance angle 
    Mtemp=[cosd(b(j,:));((sind(b(j,:))).*(1i.*q(j)));((sind(b(j,:))).*(1i./q(j)));cosd(b(j,:))];  
    Mtemp=reshape(Mtemp,2,2,numel(z)); % temporary 3D propagation matrix 
    % Since reshape works column-wise starting from the first matrix-shell
    % of a matrix and starts putting elements keeping row no. constant, 
    % Mtemp contains *q(j) term "upper" than /q(j) term in order to get
    % correct 2X2 propagation matrix.
    M(:,:,:,j)=Mtemp; 
    % Mtemp was created because reshape doesn't support changing dimension.
end
MPT=M(:,:,:,1); % initializing total propagation matrix
for k=1:numel(z)
    for j=2:N % not from j=1 since we initialized Mt from from M(:,:,:,1)
        MPT(:,:,k)=(M(:,:,k,j))*(MPT(:,:,k));
        % If..end for Total Characteristic Matrix calculation
        if (k==numel(z)) && (j==N-1)
            % propagation matrices multiplied upto only j=N-1 because analyte layer is not included to calculate reflectivity 
            MCT=inv(MPT(:,:,numel(z)));
        end
    end
end
%__________________________________________________________________________
q0=(cosd(RA)).*((sqrt(u))/np);% qj for prism at resonance angle
%__________________________________________________________________________
% REFLECTION COEFFICIENT, FIELD CALCULATION; PLOTTING
%__________________________________________________________________________
r=((MCT(1,1)+MCT(1,2).*(q(N))).*q0-(MCT(2,1)+MCT(2,2).*(q(N))))./((MCT(1,1)+MCT(1,2).*(q(N))).*q0+(MCT(2,1)+MCT(2,2).*(q(N)))); % Reflection Coefficient
Exjz=-1.*((MPT(2,1,:)).*((1+r)/q0)+(MPT(2,2,:)).*(1-r)); % x-Electric field normalized wrt. incident |Ex|.
Exjz=reshape(Exjz,1,numel(z)); % reshaped to make plottable
Ex2=(abs(Exjz)).^2; % x-Electric field modulus squared
h1=plot(z,Ex2,'.-'); % to apply legend only for this plot, mark this some variable(say h1).
%__________________________________________________________________________
% BLOCK TO PLOT VERTICAL LINES AT EVERY BOUNDARY 
% can be removed without any change to the rest of the code
%__________________________________________________________________________
hold on
y=get(gca,'ylim'); % gets current axes y limits
sumd=0;
for j=1:N-1 % N layers have N-1 boundaries
    sumd=sumd+d(j);
    h2=plot([sumd sumd],y,'m--');
    set(h2,'handlevisibility','off') % suppresses legends for all the vlines
    hold on
end
hold off
%__________________________________________________________________________
% BLOCK TO LABEL THE GRAPH
% can be removed without any change to the rest of the code
%__________________________________________________________________________
xlabel('Distance from prism interface (m)','Interpreter','latex'); 
% for normalized electric field, ylabel is as follows:
ylabel('$\displaystyle\frac{|E_x|^2}{|E^{inc}_x|^2}$','Interpreter','latex'); 
legend(h1,{strcat('Resonance Angle= ',num2str(RA),'$^\circ$')},'Location','northeast','Interpreter','latex');
%__________________________________________________________________________
toc % "tic...toc" tells runtime
%__________________________________________________________________________
% FUNCTION DECLARATION
%__________________________________________________________________________
% costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%nj varies in layers 
function bofj = bj(theta,nj,zj)
    global lambda np 
    bofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((360*nj*zj)./lambda);
    %Beta j used in Shalabney paper, %360->2*pi when cosd->cos, sind->sin
end
function qofj = qj(theta,nj)
    global u np
    qofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((sqrt(u))/nj);%qj in analyte %used cos of angle in that layer
end 