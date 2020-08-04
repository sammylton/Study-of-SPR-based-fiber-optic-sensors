clear;
N=4;
d=1e-9.*[10;40;10];%%thickness (so that Rmin<0.01, acc. to me it is making Rmin nearly vanish)of metal layer for wavelength=633nm approx. 
z=0:1e-9:2*(d(1)+d(2)+d(3));
% Block for generating separate Zj for each layer starts.
Z=zeros(N-1,numel(z));
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
