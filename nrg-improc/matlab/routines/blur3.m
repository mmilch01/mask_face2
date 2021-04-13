% Author: 	Mikhail Milchenko, mmilch@npg.wustl.edu, Washington University School of Medicine. 
% Warranty:	No express or implied warranty, or fitness of this code for specific purpose is assumed. 
% License:	You are free to use and re-distribute this code in any form, preserving the above notice.

function [B] = blur3(A,k_size)
if(rem(k_size,2)~=0) 
    k_size=k_size-1;
end;
sz=size(A);
%szx=sz(1)+2*k_size;
%szy=sz(2)+2*k_size;
%szz=sz(3)+2*k_size;
%B1=zeros(szx,szy,szz);
%B1(k_size+1:szx-k_size,k_size+1:szy-k_size,k_size+1:szz-k_size)=A;

filter=(1.0/(k_size*k_size*k_size))*double(ones(k_size,k_size,k_size));
fA=fftn(A);
ff=fftn(filter,size(A));
%B=sh(ifftn(fA.*ff),k_size/2);
B=ifftn(fA.*ff);
B=sh(B,k_size/2);

%B=zeros(sz(1),sz(2),sz(3));
%B=B1(k_size+1:szx-k_size,k_size+1:szy-k_size,k_size+1:szz-k_size);
return;
function [C]=sh(B,s)
sz=size(B);
C=double(zeros(sz));
for z=1:sz(3)
    for y=1:sz(2)
        for x=1:sz(1)
            C(x,y,z)=B(mod(x+s,sz(1))+1,mod(y+s,sz(2))+1,mod(z+s,sz(3))+1);
        end;
    end;
end;
return;
