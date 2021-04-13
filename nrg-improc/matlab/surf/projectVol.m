% Author: 	Mikhail Milchenko, mmilch@npg.wustl.edu, Washington University School of Medicine. 
% Warranty:	No express or implied warranty, or fitness of this code for specific purpose is assumed. 
% License:	You are free to use and re-distribute this code in any form, preserving the above notice.

function[Vres, Vmask]=projectVol(V,V_orig,surf,lower_surf,upper_surf,thickness,step,dir,opt)
[wx,wy,wz]=size(V_orig);
sz=size(V_orig);
%step=thickness*2;
dx=floor((wx-1)/step(1));
dy=floor((wy-1)/step(2));
%global last_tetr;

%project bounded subvolume from original volume to result volume
if(strcmp(dir,'direct')~=0) 
    Vmask=0;
    Vres=zeros(wx,wy,thickness*2+1);
    %Vres=zeros(wx,wy,wz);
    for x=1:dx
        lx=(x-1)*step(1)+1;
        rx=lx+step(1);
        for y=1:dy
          ty=(y-1)*step(2)+1;
          by=ty+step(2);
          face_b=makeFace(x,y,step,lower_surf);
          face_o=makeFace(x,y,step,surf);
          face_t=makeFace(x,y,step,upper_surf);
          
          prisms=getPrismoids(face_b,face_t,face_o);
          if(~isempty(prisms))
              el=project_prism_direct(V,prisms,step,thickness);
              Vres(lx:rx-1,ty:by-1,:)=el;
          end;
          %test
          %break;
        end;
        %test
        %break;
    end;
else  %project subvolume to bounded subvolume in original volume.
    Vres=V_orig;
    Vmask=0*Vres;
    for x=1:dx
        lx=(x-1)*step(1)+1;
        rx=lx+step(1);
        for y=1:dy
          ty=(y-1)*step(2)+1;
          by=ty+step(2);
          [ptl_l,pbl_l,ptr_l,pbr_l]=make3DPoints(x,y,step,lower_surf);
          [ptl_u,pbl_u,ptr_u,pbr_u]=make3DPoints(x,y,step,upper_surf);

          face_b=makeFace(x,y,step,lower_surf);
          face_o=makeFace(x,y,step,surf);
          face_t=makeFace(x,y,step,upper_surf);
          prisms=getPrismoids(face_b,face_t,face_o);
          
          %volume boundary.
          btl_l=min([ptl_l;pbl_l;ptr_l;pbr_l;ptl_u;pbl_u;ptr_u;pbr_u]);
          bbr_u=max([ptl_l;pbl_l;ptr_l;pbr_l;ptl_u;pbl_u;ptr_u;pbr_u]);
          btl_l=atb(floor(btl_l),[wx,wy,wz],3);
          bbr_u=atb(ceil(bbr_u),[wx,wy,wz],3);
          if(~isempty(prisms))
              orig=[lx,ty];
%[Vres, Vmask] = project_prism_reverse1(V,Vres,Vmask,prisms,[lx,ty],step,thickness);              
%here goes the code for project_prism_reverse1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Vres=Vr;
%            Vmask=Vm;
            tetra_old(1:3)=prism_partition(makeRegularPrismoid(prisms(1).ind,orig,step,thickness));
            tetra_new(1:3)=prism_partition(prisms(1));

            pmin1=min(min(prisms(1).bt,prisms(1).tp));
            pmax1=max(max(prisms(1).bt,prisms(1).tp));

            if(size(prisms,2)>1)
                tetra_old(4:6)=prism_partition(makeRegularPrismoid(prisms(2).ind,orig,step,thickness));
                tetra_new(4:6)=prism_partition(prisms(2));
                pmin2=min(min(prisms(2).bt,prisms(2).tp));
                pmax2=max(max(prisms(2).bt,prisms(2).tp));   
            else
                pmin2=pmin1;
                pmax2=pmax1;
            end;
            pmin=round(min(sz,max([1,1,1],min(pmin1,pmin2))));
            pmax=round(max([1,1,1],min(sz,max(pmax1,pmax2))));

            %main cycle.
            t=1;
            for zz=pmin(3):pmax(3)
                for yy=pmin(2):pmax(2)
                    for xx=pmin(1):pmax(1)
                        pt=[xx,yy,zz];%-[1,1,1];                 
                        t=get_containing_tetrahedron(pt,tetra_new,t);
                        if(t<1)
%                            if(Vmask(xx,yy,zz)==0); Vmask(xx,yy,zz)=5; end;
                            t=1;
                            continue;
                        end;
                        x_old=round(coord_transform(pt,tetra_new(t),tetra_old(t)));
                        %test reverse transform            
                        x_old1=atb_l(x_old,[1,1,1],sz,3);
                        if(norm(x_old1-x_old)>0) %should not happen
                            continue;
                        end;
                        val=V(x_old1(1),x_old1(2),x_old1(3));        
                        Vres(xx,yy,zz)=val;
                        Vmask(xx,yy,zz)=20+2*t;
                    end;
                end;
            end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end of project_prism_reverse1                            
          end;                               
          %test
          %break;
        end;
        %test
        %break;
    end;    
end;
    
function [face] = makeFace(x,y,step,surf)
lx=(x-1)*step(1)+1;
rx=lx+step(1);
ty=(y-1)*step(2)+1;
by=ty+step(2);

vert(1,:)=surf(1:3,x,y)';
vert(2,:)=surf(1:3,x+1,y)';
vert(3,:)=surf(1:3,x+1,y+1)';
vert(4,:)=surf(1:3,x,y+1)';


face.vert=vert;
deg=ones(4,1);
for(i=1:4) 
   if(vert(i,3)~=1) deg(i)=0;end;
end;
face.deg=deg;
face.ndeg=sum(deg);

function [pr] = makePrismoid(face_b,face_t,ind)
for i=1:3
    bt(i,:)=face_b.vert(ind(i),:);
    tp(i,:)=face_t.vert(ind(i),:);
end;
pr.bt=bt;
pr.tp=tp;
pr.ind=ind;

function [pr] = makeRegularPrismoid(ind,origin,step,thickness)
lx=origin(1);ty=origin(2);rx=origin(1)+step(1)-1;by=origin(2)+step(2)-1;
nzmax=thickness*2+1;

ptl(1,:)=[lx,ty,1]; ptu(1,:)=[lx,ty,nzmax];
ptl(2,:)=[rx,ty,1]; ptu(2,:)=[rx,ty,nzmax];
ptl(3,:)=[rx,by,1]; ptu(3,:)=[rx,by,nzmax];
ptl(4,:)=[lx,by,1]; ptu(4,:)=[lx,by,nzmax];

for i=1:3
    pb(i,:)=ptl(ind(i),:);
    pt(i,:)=ptu(ind(i),:);
end;
pr.bt=pb;
pr.tp=pt;
pr.ind=ind;
    
    
function [prisms] = getPrismoids(face_b,face_t,face_orig)
if(face_orig.ndeg>1) prisms=[]; return; end;
deg=face_orig.deg;
if(face_orig.ndeg==1) %one vertex is degenerate, exclude containing prismoid
    if(deg(1)==1) ind=[2,3,4];
    elseif(deg(2)==1) ind=[3,4,1];
    elseif(deg(3)==1) ind=[4,1,2];
    else ind=[1,2,3];
    end;
    prisms(1)=makePrismoid(face_b,face_t,ind);
else %regular case.
    prisms(1)=makePrismoid(face_b,face_t,[1,2,3]);
    prisms(2)=makePrismoid(face_b,face_t,[1,4,3]);
end;

function [tl, bl, tr, br]=make3DPoints(x,y,step,surf)
lx=(x-1)*step(1)+1;
rx=lx+step(1);
ty=(y-1)*step(2)+1;
by=ty+step(2);
tl=[lx,ty,surf(3,x,y)];
bl=[lx,by,surf(3,x,y+1)];
tr=[rx,ty,surf(3,x+1,y)];
br=[rx,by,surf(3,x+1,y+1)];

function [Vres,Vmask]=project_prism_reverse1(V,Vr,Vm,prisms,orig,step,thickness)
sz=size(Vr);
Vres=Vr;
Vmask=Vm;
tetra_old(1:3)=prism_partition(makeRegularPrismoid(prisms(1).ind,orig,step,thickness));
tetra_new(1:3)=prism_partition(prisms(1));

pmin1=min(min(prisms(1).bt,prisms(1).tp));
pmax1=max(max(prisms(1).bt,prisms(1).tp));

if(size(prisms,2)>1)
    tetra_old(4:6)=prism_partition(makeRegularPrismoid(prisms(2).ind,orig,step,thickness));
    tetra_new(4:6)=prism_partition(prisms(2));
    pmin2=min(min(prisms(2).bt,prisms(2).tp));
    pmax2=max(max(prisms(2).bt,prisms(2).tp));   
else
    pmin2=pmin1;
    pmax2=pmax1;
end;
pmin=round(min(sz,max([1,1,1],min(pmin1,pmin2))));
pmax=round(max([1,1,1],min(sz,max(pmax1,pmax2))));

%main cycle.
t=1;
for z=pmin(3):pmax(3)
    for y=pmin(2):pmax(2)
        for x=pmin(1):pmax(1)
            pt=[x,y,z];%-[1,1,1];                 
            t=get_containing_tetrahedron(pt,tetra_new,t);
            if(t<1)
                t=1;
                if(Vmask(x,y,z)==0); Vmask(x,y,z)=5; end;
                continue;
            end;
            
            x_old=round(coord_transform(pt,tetra_new(t),tetra_old(t)));

            %test reverse transform            
            x_old1=atb_l(x_old,[1,1,1],sz,3);
            if(norm(x_old1-x_old)>0) %should not happen
                continue;
            end;
            val=V(x_old1(1),x_old1(2),x_old1(3));        
            Vres(x,y,z)=val;
            Vmask(x,y,z)=20+2*t;
        end;
    end;
end;


function [Vres,Vmask]=project_prism_reverse(V,V_orig,prisms,orig,btl_l,bbr_u,step,thickness)
sz=size(V);
nzmax=thickness*2+1;
tetra_old(1:3)=prism_partition(makeRegularPrismoid(prisms(1).ind,orig,step,thickness));
tetra_new(1:3)=prism_partition(prisms(1));

if(size(prisms,2)>1)
    tetra_old(4:6)=prism_partition(makeRegularPrismoid(prisms(2).ind,orig,step,thickness));
    tetra_new(4:6)=prism_partition(prisms(2));
end;

mx=bbr_u-btl_l+[1,1,1];

Vres=V_orig(btl_l(1):bbr_u(1),btl_l(2):bbr_u(2),btl_l(3):bbr_u(3));
Vmask=0*Vres;

bound_low=[orig(1),orig(2),1];
bound_high=bound_low+[step,step,nzmax];
%main cycle.
t=1;
for z=1:mx(3)
    for y=1:mx(2)
        for x=1:mx(1)
            pt=[x,y,z]+btl_l-[1,1,1];                 
            t=get_containing_tetrahedron(pt,tetra_new,t);
            if(t<1)
                t=1;
                Vmask(x,y,z)=10;
                continue;
            end;
            
            x_old=round(coord_transform(pt,tetra_new(t),tetra_old(t)));

            %test reverse transform            
            x_old1=atb_l(x_old,bound_low,bound_high,3);
            if(norm(x_old1-x_old)>0) %should not happen
                continue;
            end;
            val=V(x_old1(1),x_old1(2),x_old1(3));        
            Vres(x,y,z)=val;
            Vmask(x,y,z)=20+2*t;
        end;
    end;
end;

function [Vres]=project_prism_direct(V,prisms,step,thickness)
      
sz=size(V);
tetra_old(1:3)=prism_partition(prisms(1));
tetra_new(1:3)=prism_partition(makeRegularPrismoid(prisms(1).ind,[1,1],step,thickness));

if(size(prisms,2)>1)
    tetra_old(4:6)=prism_partition(prisms(2));
    tetra_new(4:6)=prism_partition(makeRegularPrismoid(prisms(2).ind,[1,1],step,thickness));
end;

nzmax=thickness*2+1;
Vres=zeros(step(1),step(2),nzmax);

%main cycle.
t=1;
for z=1:nzmax
    for y=1:step(2)
        for x=1:step(1)
            t=get_containing_tetrahedron([x,y,z],tetra_new,t);
            if(t<1)
                t=1;
                continue;
            end;
            
            x_old=round(coord_transform([x,y,z],tetra_new(t),tetra_old(t)));
            x_old1=atb(x_old,sz,3);
            val=V(x_old1(1),x_old1(2),x_old1(3));        
            Vres(x,y,z)=20+2*t+0.01;
            Vres(x,y,z)=val;
        end;
    end;
end;

function[t]=get_containing_tetrahedron(x,tetra,last_tetr)
%global last_tetr;
t=-1;
ntetr=size(tetra,2);
if (last_tetr>ntetr) last_tetr=1; end;
for i=1:ntetr
    k=i;
    if(i==1) k=last_tetr; end;
    if(i==last_tetr) k=1; end;      
    if(tetra(k).deg~=0) res=0; continue; end;
    res=1;
    for j=1:4
        if(x(1)*tetra(k).faces(1,j)+x(2)*tetra(k).faces(2,j)+x(3)*tetra(k).faces(3,j)+tetra(k).faces(4,j)<0)
        %   if(ev_plane_eq(x,tetr.faces(:,i))<dscr)
            res=0; break;
        end;
    end;    
    if(res>0) t=k; return; end;
%   if(is_inside_tetrahedron(x,tetra(i))>0) t=i;return;end;
end;

function[res]=is_inside_tetrahedron(x,tetr)
dscr=-1e-10;
if(tetr.deg~=0) res=0; return; end;
res=1;
for i=1:4
    if(x(1)*tetr.faces(1,i)+x(2)*tetr.faces(2,i)+x(3)*tetr.faces(3,i)+tetr.faces(4,i)<0)
%   if(ev_plane_eq(x,tetr.faces(:,i))<dscr)
        res=0; return;
    end;
end;

function[tetra]=prism_partition(prism)
tetra(1)=tetrahedron(prism.bt(3,:),prism.bt(1,:),prism.bt(2,:),prism.tp(1,:));
tetra(2)=tetrahedron(prism.tp(2,:),prism.tp(1,:),prism.tp(3,:),prism.bt(3,:));
tetra(3)=tetrahedron(prism.tp(1,:),prism.tp(2,:),prism.bt(2,:),prism.bt(3,:));


function[coefs]=plane_eq(x1,x2,x3)
cx=(x1(2)-x2(2))*(x3(3)-x2(3))-(x1(3)-x2(3))*(x3(2)-x2(2));
cy=(x1(3)-x2(3))*(x3(1)-x2(1))-(x1(1)-x2(1))*(x3(3)-x2(3));
cz=(x1(1)-x2(1))*(x3(2)-x2(2))-(x1(2)-x2(2))*(x3(1)-x2(1));
c0=(-x2(2)*x3(1)+x2(1)*x3(2))*  x1(3)+...
   ( x2(3)*x3(1)-x2(1)*x3(3))*  x1(2)+...
   (-x2(3)*x3(2)+x2(2)*x3(3))* x1(1);
coefs=-[cx,cy,cz,c0];
mx=max(coefs);
if(mx>0) coefs=coefs/mx; end;


function[res]=tetrahedron(x1,x2,x3,x4)
coefs=zeros(4,4);
base=plane_eq(x1,x2,x3);
left=plane_eq(x2,x4,x3);
right=plane_eq(x4,x2,x1);
rear=plane_eq(x1,x3,x4);
coefs(:,1)=base;
coefs(:,2)=left;
coefs(:,3)=right;
coefs(:,4)=rear;

%equation sign checks.
sgn=ones(1,4);
if(x4(1)*base(1)+x4(2)*base(2)+x4(3)*base(3)+base(4)<0)
%if (ev_plane_eq(x4,base)<0)
    sgn(1)=-1;
end;
if(x1(1)*left(1)+x1(2)*left(2)+x1(3)*left(3)+left(4)<0)
%if (ev_plane_eq(x1,left)<0)
    sgn(2)=-1;
end;
if(x3(1)*right(1)+x3(2)*right(2)+x3(3)*right(3)+right(4)<0)
%if (ev_plane_eq(x3,right)<0)
    sgn(3)=-1;
end;
if(x2(1)*rear(1)+x2(2)*rear(2)+x2(3)*rear(3)+rear(4)<0)
%if(ev_plane_eq(x2,rear)<0)
    sgn(4)=-1;
end;
for i=1:4
    coefs(:,i)=coefs(:,i)*sgn(i);
end;

res.sgn=sgn;
res.faces=coefs;
res.origin=x2;
v1=x1-x2;v2=x3-x2;v3=x4-x2;
n1=norm(v1);n2=norm(v2);n3=norm(v3);

%if(n1>0)v1=v1/n1;end;
%if(n2>0)v2=v2/n2;end;
% if(n3>0)v3-v3/n3;end;
res.v1=v1;res.v2=v2;res.v3=v3;
res.n1=n1;res.n2=n2;res.n3=n3;

if(n1>0 && n2>0 && n3>0)
    res.deg=0;
else res.deg=1;
end;
%useful, dot parallel products.
%dp1=(left(1)*v1(1)+left(2)*v1(2)+left(3)*v1(3));
%dp2=(right(1)*v2(1)+right(2)*v2(2)+right(3)*v2(3));
%dp3=(base(1)*v3(1)+base(2)*v3(2)+base(3)*v3(3));
dp1=(coefs(1,2)*v1(1)+coefs(2,2)*v1(2)+coefs(3,2)*v1(3));
dp2=(coefs(1,3)*v2(1)+coefs(2,3)*v2(2)+coefs(3,3)*v2(3));
dp3=(coefs(1,1)*v3(1)+coefs(2,1)*v3(2)+coefs(3,1)*v3(3));

res.dp=[dp1,dp2,dp3];
n=zeros(3,3);
n(:,1)=n1;
n(:,2)=n2;
n(:,3)=n3;
res.n=n;

%vx=[v1(1),v2(1),v3(1)];
%vy=[v1(2),v2(2),v3(2)];
%vz=[v1(3),v2(3),v3(3)];
%cx=[coefs(1,2),coefs(1,3),coefs(1,1)];
%cy=[coefs(2,2),coefs(2,3),coefs(2,1)];
%cz=[coefs(3,2),coefs(3,3),coefs(3,1)];
%res.vx=vx;res.vy=vy;res.vz=vz;
%res.cx=cx;res.cy=cy;res.cz=cz;

function[res]=ev_plane_eq(x,eq)
res=eq(1)*x(1)+eq(2)*x(2)+eq(3)*x(3)+eq(4);

function[xr]=atb_l(x,sz_l,sz_u,nDim)
xr=x;
for i=1:nDim
    if(x(i)<sz_l(i)) xr(i)=sz_l(i); end;
    if(x(i)>sz_u(i)) xr(i)=sz_u(i);end;
end;

function[xr]=atb(x,sz,nDim)
xr=x;
for i=1:nDim
    if(x(i)<1) xr(i)=1; end;
    if(x(i)>sz(i)) xr(i)=sz(i);end;
end;

%project point x on axis in tetrahedral coord system.
function[r]=project_on_axis(xs,num,tetr)
x=xs(1);y=xs(2);z=xs(3);
if(num==1)
    n=tetr.n1;
    dp=tetr.dp(1);
    if(tetr.dp(1)<=0) return; end;
    vx=tetr.v1(1);
    vy=tetr.v1(2);
    vz=tetr.v1(3);
    cx=tetr.faces(1,2);
    cy=tetr.faces(2,2);
    cz=tetr.faces(3,2);
elseif(num==2)
    n=tetr.n2;
    dp=tetr.dp(2);
    if(tetr.dp(2)<=0) return; end;
    vx=tetr.v2(1);
    vy=tetr.v2(2);
    vz=tetr.v2(3);
    cx=tetr.faces(1,3);
    cy=tetr.faces(2,3);
    cz=tetr.faces(3,3);
else %num==3
    n=tetr.n3;
    dp=tetr.dp(3);
    if(tetr.dp(3)<=0) return; end;
    vx=tetr.v3(1);
    vy=tetr.v3(2);
    vz=tetr.v3(3);
    cx=tetr.faces(1,1);
    cy=tetr.faces(2,1);
    cz=tetr.faces(3,1);
end;
x0=tetr.origin(1);y0=tetr.origin(2);z0=tetr.origin(3);
xp=((vx*z+x0*vz-vx*z0)*cz+(vx*y+vy*x0-vx*y0)*cy+vx*cx*x)/dp;
yp=((z*vy+vz*y0-vy*z0)*cz+(x*vy-vy*x0+vx*y0)*cx+cy*y*vy)/dp;
zp=((-vz*y0+vz*y+vy*z0)*cy+(-x0*vz+vz*x+vx*z0)*cx+vz*cz*z)/dp;
r=norm([xp-x0,yp-y0,zp-z0])/n;
return;

function[x_dest]=coord_transform(xs,tetr,td)
x=xs(1);y=xs(2);z=xs(3);
x0=tetr.origin(1);y0=tetr.origin(2);z0=tetr.origin(3);
%r1
n=tetr.n1;
dp=tetr.dp(1);
if(tetr.dp(1)<=0) 
    r1=0; 
else    
    vx=tetr.v1(1);
    vy=tetr.v1(2);
    vz=tetr.v1(3);
    cx=tetr.faces(1,2);
    cy=tetr.faces(2,2);
    cz=tetr.faces(3,2);
    xp=((vx*z+x0*vz-vx*z0)*cz+(vx*y+vy*x0-vx*y0)*cy+vx*cx*x)/dp;
    yp=((z*vy+vz*y0-vy*z0)*cz+(x*vy-vy*x0+vx*y0)*cx+cy*y*vy)/dp;
    zp=((-vz*y0+vz*y+vy*z0)*cy+(-x0*vz+vz*x+vx*z0)*cx+vz*cz*z)/dp;
    r1=norm([xp-x0,yp-y0,zp-z0])/n;            
%    r1=project_on_axis(x_src,1,ts);
end;
n=tetr.n2;
dp=tetr.dp(2);
if(tetr.dp(2)<=0) 
    r2=0;
else    
    vx=tetr.v2(1);
    vy=tetr.v2(2);
    vz=tetr.v2(3);
    cx=tetr.faces(1,3);
    cy=tetr.faces(2,3);
    cz=tetr.faces(3,3);
    xp=((vx*z+x0*vz-vx*z0)*cz+(vx*y+vy*x0-vx*y0)*cy+vx*cx*x)/dp;
    yp=((z*vy+vz*y0-vy*z0)*cz+(x*vy-vy*x0+vx*y0)*cx+cy*y*vy)/dp;
    zp=((-vz*y0+vz*y+vy*z0)*cy+(-x0*vz+vz*x+vx*z0)*cx+vz*cz*z)/dp;
    r2=norm([xp-x0,yp-y0,zp-z0])/n;            
%   r2=project_on_axis(x_src,2,ts);
end;
n=tetr.n3;
dp=tetr.dp(3);
if(tetr.dp(3)<=0)
    r3=0;
else
    vx=tetr.v3(1);
    vy=tetr.v3(2);
    vz=tetr.v3(3);
    cx=tetr.faces(1,1);
    cy=tetr.faces(2,1);
    cz=tetr.faces(3,1);
    xp=((vx*z+x0*vz-vx*z0)*cz+(vx*y+vy*x0-vx*y0)*cy+vx*cx*x)/dp;
    yp=((z*vy+vz*y0-vy*z0)*cz+(x*vy-vy*x0+vx*y0)*cx+cy*y*vy)/dp;
    zp=((-vz*y0+vz*y+vy*z0)*cy+(-x0*vz+vz*x+vx*z0)*cx+vz*cz*z)/dp;
    r3=norm([xp-x0,yp-y0,zp-z0])/n;                
%    r3=project_on_axis(x_src,3,ts);
end;
x_dest=td.origin+td.v1*r1+td.v2*r2+td.v3*r3;
