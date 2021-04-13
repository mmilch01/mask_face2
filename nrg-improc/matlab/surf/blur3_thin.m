% Send comments to Mikhail Milchenko, mmilch@wustl.edu

function [B] = blur3_thin(A,thickness)

sz=size(A);

fz=max(2,ceil(thickness/2));
fx=ceil(sz(1)*.1); 
fy=ceil(sz(2)*.1);
A1=pad(A,fx,fy,fz);

f=(1/((fx*2+1)*(fy*2+1)*(fz*2+1)))*ones(fx*2+1,fy*2+1,fz*2+1);
B=convn(A1,f,'valid');

minxy=max(2,ceil(thickness/2));
minz=max(1,ceil(thickness*.05));
rx=fx-minxy; ry=fy-minxy; rz=fz-minz;
B1=B;
%adjust the lower layer.
for z=1:rz
    wz=round(minz+rz*sin((z-1)*pi/(2*rz)));
    wx=round(minxy+rx*sin((z-1)*pi/(2*rz)));
    wy=round(minxy+ry*sin((z-1)*pi/(2*rz)));
    f=(1/((wx*2+1)*(wy*2+1)*(wz*2+1)))*ones(wx*2+1,wy*2+1,wz*2+1);
    A1=pad(A,wx,wy,wz);
    B0=convn(A1,f,'valid');
    B1(:,:,z)=B0(:,:,z);
end    
B=B1;

% for x=1:sz(1)
%     for y=1:sz(2)
%         lowval=A(x,y,1);
%         highval=A(x,y,dz);
%         for z=1:dz
%           B(x,y,z)=lowval+(highval-lowval)*(z/dz);
%         end;
%     end;
% end;
function[B]=pad(A,fx,fy,fz)
sz=size(A);
B=zeros(sz(1)+2*fx,sz(2)+2*fy,sz(3)+2*fz);
for z=1-fz:sz(3)+fz
    for y=1-fy:sz(2)+fy
        for x=1-fx:sz(1)+fx
            x1=atb_l([x,y,z],[1,1,1],[sz(1),sz(2),sz(3)],3);
            B(x+fx,y+fy,z+fz)=A(x1(1),x1(2),x1(3));
        end;
    end;
end;
            
function[xr]=atb_l(x,sz_l,sz_u,nDim)
xr=x;
for i=1:nDim
    if(x(i)<sz_l(i)) xr(i)=sz_l(i); end;
    if(x(i)>sz_u(i)) xr(i)=sz_u(i);end;
end;
            