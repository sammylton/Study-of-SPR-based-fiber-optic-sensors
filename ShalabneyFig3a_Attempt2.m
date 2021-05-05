clear;
lambda=633e-9;
nm=0.1325+4.0203i;%real part is normal refractive index and complex part
%is the decay term
%nm=0.1325;
erm=nm^2;
np=1.732;
na=1.33;
u=pi*4e-7;
dm=43e-9; 
theta=linspace(0,90,1000000);
costheta=sqrt(1-((sind(theta)).^2).*((np^2)/(nm^2)));
b=(costheta).*((360*nm*dm)/lambda);
q0=(cosd(theta)).*((sqrt(u))/np);
qn=(sqrt(1-((sind(theta)).^2).*((np^2)/(na^2)))).*((sqrt(u))/na);
q=costheta.*(sqrt(u/erm));
%M=[cos(b),-((sin(b)).*(1i./q));-((sin(b)).*(1i.*q)),cos(b)];
M11=cosd(b);
M12=-((sind(b)).*(1i./q));
M21=-((sind(b)).*(1i.*q));
M22=cosd(b);
%r=((M(1,1)+(M(1,2)).*qn).*q0-(M(2,1)+(M(2,2)).*qn))./((M(1,1)+(M(1,2)).*qn).*q0+(M(2,1)+(M(2,2)).*qn));
r=((M11+M12.*qn).*q0-(M21+M22.*qn))./((M11+M12.*qn).*q0+(M21+M22.*qn));
R=(abs(r)).^2;
plot(theta,R);
%Below is the block to find extrema
indexmin = find(min(R) == R); 
thetamin = theta(indexmin); 
Rmin = R(indexmin);
indexmax = find(max(R) == R);
thetamax = theta(indexmax); 
Rmax = R(indexmax); 
%Above is the black to find extrema