% For normalised electric field. In ShalabneyFig4a4, I used different matrix manipulations.
% This code changes those manipulations so that a foundational code is made for multiple layer case 
% so that matrix inversion/multiplication is a straight-forward step.
% And this code could contain resonance angle finding block
% One can also use the metal thickness optimisation code block here from FindMinThickness.m
clear all;
tic
% Shalabney Fig6b
u=pi*4e-7;%magnetic permeability of the assumed non-magnetic substances
lambda=633e-9;%wavelength of incident light in vacuum at which nm is as follows
np=1.732;
na=1.33;
N=3;
n=[0.1325+4.0203*1i;3.8354+0.0245*1i;na];%refractive indices of layers including that of analyte at wavelength=633nmop
for o=1:5
    if o==1, dloop=0; end
    if o==2, dloop=4; end
    if o==3, dloop=8; end
    if o==4, dloop=10.5; end
    if o==5, dloop=12; end
    
    d=1e-9.*[43,dloop];
    
    %_________________________________________________________________________________________________________________________________________________________
    %block to find resonance theta for given values starts
    %_________________________________________________________________________________________________________________________________________________________
    theta=0:1e-3:90;%angle in prism(angle of incidence)
    Mtot=zeros(2,2,numel(theta));
    
    %initialization loop
    for l=1:numel(theta)
        %Preallocationg M in this loop for the next loop
        Mtot(:,:,l)=eye(2,2);
    end
    
    for j=1:N-1
        %Mk=[cosd(bj(theta,n(j),d(j)));-((sind(bj(theta,n(j),d(j)))).*(1i./qj(theta,n(j))));-((sind(bj(theta,n(j),d(j)))).*(1i.*qj(theta,n(j))));cosd(bj(theta,n(j),d(j)))];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
        Mk=[cosd(bj(theta,n(j),d(j)));((sind(bj(theta,n(j),d(j)))).*(-1i.*qj(theta,n(j))));((sind(bj(theta,n(j),d(j)))).*(-1i./qj(theta,n(j))));cosd(bj(theta,n(j),d(j)))];%Inverse(hence 1,2 and 2,1 elements are negated) of complete propagation matrix in metal
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
    rdash=((M11+(M12).*(qj(theta,na))).*(qj(theta,np))-(M21+(M22).*(qj(theta,na))))./((M11+(M12).*(qj(theta,na))).*(qj(theta,np))+(M21+(M22).*(qj(theta,na))));%reflection coefficient
    R=(abs(rdash)).^2;%reflectivity
    %Below is the block to find extrema
    indexmin = find(min(R) == R);
    resonancetheta = theta(indexmin);
    %_________________________________________________________________________________________________________________________________________________________
    %block to find resonance theta for given values ends
    %_________________________________________________________________________________________________________________________________________________________
    
    z=0:0.5e-9:270e-9;%for x-axis of required graph
    % 0.5e-9 because we want to get 10.5e-9 in the data too for precise boundary
    z1=z;
    z1(z>=d(1))=d(1);%> or >= makes no difference
    
    z2=z-d(1);
    z2(z<=d(1))=0;%> or >= makes no difference
    z2(z>=d(1)+d(2))=d(2);
    
    z3=z-(d(1)+d(2));
    z3(z<=d(1)+d(2))=0;
    
    q0=(cosd(resonancetheta)).*((sqrt(u))/np);
    
    costheta1=sqrt(1-((sind(resonancetheta)).^2).*((np^2)/((n(1))^2)));%angle in metal
    costheta2=sqrt(1-((sind(resonancetheta)).^2).*((np^2)/((n(2))^2)));%angle in dielectric
    costheta3=sqrt(1-((sind(resonancetheta)).^2).*((np^2)/((n(3))^2)));%angle in analyte
    b1=(costheta1).*((360*(n(1))/lambda).*z1);
    b2=(costheta2).*((360*(n(2))/lambda).*z2);
    b3=(costheta3).*((360*(n(3))/lambda).*z3);
    q1=costheta1.*(sqrt(u)/n(1));
    q2=costheta2.*(sqrt(u)/n(2));
    q3=costheta3.*(sqrt(u)/n(3));
    % Matrix for propagation in metal
    M1=[cosd(b1);(sind(b1).*(1i.*q1));(sind(b1).*(1i./q1));cosd(b1)];
    M1=reshape(M1,2,2,numel(z));
    % Matrix for propagation in dielectric
    M2=[cosd(b2);(sind(b2).*(1i.*q2));(sind(b2).*(1i./q2));cosd(b2)];
    M2=reshape(M2,2,2,numel(z));
    % Matrix for propagation in analyte
    M3=[cosd(b3);(sind(b3).*(1i.*q3));(sind(b3).*(1i./q3));cosd(b3)];
    M3=reshape(M3,2,2,numel(z));
    Mtemp=zeros(2,2,numel(z));
    for k=1:numel(z)
        Mtemp(:,:,k)=(M2(:,:,k))*(M1(:,:,k));
    end
    % Matrix for calculation of r
    Mr=inv(Mtemp(:,:,numel(z)));
    
    M=zeros(2,2,numel(z));
    for k=1:numel(z)
        M(:,:,k)=(M3(:,:,k))*(M2(:,:,k))*(M1(:,:,k));
    end
    r=((Mr(1,1)+Mr(1,2).*q3).*q0-(Mr(2,1)+Mr(2,2).*q3))./((Mr(1,1)+Mr(1,2).*q3).*q0+(Mr(2,1)+Mr(2,2).*q3));
    Exjz=-1.*((M(2,1,:)).*((1+r)/q0)+(M(2,2,:)).*(1-r));
    Exjz=reshape(Exjz,1,numel(z));
    modEx2=(abs(Exjz)).^2;
    znm=z.*(1e9);%for x-axis to be nm scale
    plot(znm,modEx2,'.-');%find a general method to create a vertical line
    hold on;
    if o==4
        %PLOT VLINES AT EVERY BOUNDARY
        countofvlines=0;
        sumd=0;
        for j=1:N-1
            sumd=sumd+d(j)*1e9;
            plot([sumd sumd],[0 max(modEx2)],'m--')
            hold on;
            countofvlines=countofvlines+1;
        end
        hold on;
    end
end
toc 
function bofj = bj(theta,nj,dj)
% BJ cosoftheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%cos of angle in metal
    global lambda np 
    bofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((360*nj*dj)./lambda);%Beta j used in Shalabney paper %360->2*pi when cosd->cos, sind->sin
end

function qofj = qj(theta,nj)
% QJ cosoftheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)));%cos of angle in metal
    global u np
    qofj=(sqrt(1-((sind(theta)).^2).*((np^2)/(nj^2)))).*((sqrt(u))/nj);%qj in analyte %used cos of angle in that layer
end 