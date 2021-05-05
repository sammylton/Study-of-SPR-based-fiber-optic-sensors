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
n=1000;
for k=1:n
    theta=(pi/2)*k/n;
    costheta=sqrt(1-((sin(theta))^2)*((np^2)/(nm^2)));
    b=(costheta)*((2*pi*nm*dm)/lambda);
    q0=(cos(theta))*((sqrt(u))/np);
    qn=(sqrt(1-((sin(theta))^2)*((np^2)/(na^2))))*((sqrt(u))/na);
    q=costheta*(sqrt(u/erm));
    M=[cos(b),-((sin(b))*1i)/q;-((sin(b))*1i*q),cos(b)];
    r=((M(1,1)+(M(1,2))*qn)*q0-(M(2,1)+(M(2,2))*qn))/((M(1,1)+(M(1,2))*qn)*q0+(M(2,1)+(M(2,2))*qn));
    R=(abs(r))^2;
    plot(theta*180/pi,R,'-x');
    hold on;
end
hold off;%so that after changing value of na, running program won't make graph in the same figure window
%plot(theta,r,'b');
%indexmin = find(min(R) == R); 
%thetamin = theta(indexmin); 
%Rmin = R(indexmin);
%indexmax = find(max(R) == R);
%thetamax = theta(indexmax); 
%Rmax = R(indexmax); 
