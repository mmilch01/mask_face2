function[Vmask]=maskVol(V_orig,lower_surf,orig_surf,upper_surf,step)
[wx,wy,wz]=size(V_orig);
sz=size(V_orig);
%step=thickness*2;
dx=floor((wx-1)/step(1));
dy=floor((wy-1)/step(2));

Vres=V_orig;
Vmask=0*Vres;
for x=1:dx
    for y=1:dy
         %volume boundary.
         [ptl_l,pbl_l,ptr_l,pbr_l]=make3DPoints(x,y,step,lower_surf);
         [ptl_u,pbl_u,ptr_u,pbr_u]=make3DPoints(x,y,step,upper_surf);
         btl_l=min([ptl_l;pbl_l;ptr_l;pbr_l;ptl_u;pbl_u;ptr_u;pbr_u]);
         bbr_u=max([ptl_l;pbl_l;ptr_l;pbr_l;ptl_u;pbl_u;ptr_u;pbr_u]);
         btl_l=atb(floor(btl_l),[wx,wy,wz],3);
         bbr_u=atb(ceil(bbr_u),[wx,wy,wz],3);
         
         %determine faces and prisms.
         face_b=makeFace(x,y,step,lower_surf);
         face_o=makeFace(x,y,step,orig_surf);
         face_t=makeFace(x,y,step,upper_surf);
         prisms=getPrismoids(face_b,face_t,face_o);
         
         if(~isempty(prisms))                          
             %el_mask=mask_prisms(V_orig,prisms,bbr_u,btl_l);
             %Vmask(btl_l(1):bbr_u(1),btl_l(2):bbr_u(2),btl_l(3):bbr_u(3))=el_mask;             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   start mask_prisms
            tetra_new(1:3)=prism_partition(prisms(1));
            pmin1=min(min(prisms(1).bt,prisms(1).tp));
            pmax1=max(max(prisms(1).bt,prisms(1).tp));
            if(size(prisms,2)>1) 
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
            for zz=pmin(3):pmax(3)
                for yy=pmin(2):pmax(2)
                    for xx=pmin(1):pmax(1)
                        pt=[xx,yy,zz];
                        t=get_containing_tetrahedron(pt,tetra_new);
                        if(t>0)
                            Vmask(xx,yy,zz)=1;%20+t*2;
                        end;
                    end;
                end;
            end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   end mask_prisms
         end;
         %testing
         %break;
    end;
    %testing
    %break;
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

function[Vmask]=mask_prisms(V_orig,prisms,bbr_u,btl_l)
tetra_new(1:3)=prism_partition(prisms(1));
if(size(prisms,2)>1) 
    tetra_new(4:6)=prism_partition(prisms(2));
end;
mx=bbr_u-btl_l+[1,1,1];
Vmask=0*(V_orig(btl_l(1):bbr_u(1),btl_l(2):bbr_u(2),btl_l(3):bbr_u(3)));
%main cycle.
for z=1:mx(3)
    for y=1:mx(2)
        for x=1:mx(1)
            pt=[x,y,z]+btl_l-[1,1,1];                 
            t=get_containing_tetrahedron(pt,tetra_new);
            if(t>0)
                Vmask(x,y,z)=20+t*2;
            else
                Vmask(x,y,z)=0;
            end;            
        end;
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

function[tetra]=prism_partition(prism)
tetra(1)=tetrahedron(prism.bt(3,:),prism.bt(1,:),prism.bt(2,:),prism.tp(1,:));
tetra(2)=tetrahedron(prism.tp(2,:),prism.tp(1,:),prism.tp(3,:),prism.bt(3,:));
tetra(3)=tetrahedron(prism.tp(1,:),prism.tp(2,:),prism.bt(2,:),prism.bt(3,:));


function [pr] = makePrismoid(face_b,face_t,ind)
for i=1:3
    bt(i,:)=face_b.vert(ind(i),:);
    tp(i,:)=face_t.vert(ind(i),:);
end;
pr.bt=bt;
pr.tp=tp;
pr.ind=ind;

function[t]=get_containing_tetrahedron(x,tetra)
t=-1;
nTetr=size(tetra,2);
for i=1:nTetr
    if(is_inside_tetrahedron(x,tetra(i))>0) t=i;return;end;
end;

function[coefs]=plane_eq(x1,x2,x3)
cx=(x1(2)-x2(2))*(x3(3)-x2(3))-(x1(3)-x2(3))*(x3(2)-x2(2));
cy=(x1(3)-x2(3))*(x3(1)-x2(1))-(x1(1)-x2(1))*(x3(3)-x2(3));
cz=(x1(1)-x2(1))*(x3(2)-x2(2))-(x1(2)-x2(2))*(x3(1)-x2(1));
c0=(-x2(2)*x3(1)+x2(1)*x3(2))*x1(3)+(x2(3)*x3(1)-x2(1)*x3(3))*x1(2)+...
    (-x2(3)*x3(2)+x2(2)*x3(3))*x1(1);
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
if (ev_plane_eq(x4,base)<0)
    sgn(1)=-1;
end;
if (ev_plane_eq(x1,left)<0)
    sgn(2)=-1;
end;
if (ev_plane_eq(x3,right)<0)
    sgn(3)=-1;
end;
if(ev_plane_eq(x2,rear)<0)
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

res.v1=v1;res.v2=v2;res.v3=v3;
res.n1=n1;res.n2=n2;res.n3=n3;

if(n1>0 && n2>0 && n3>0)
    res.deg=0;
else res.deg=1;
end;

dp1=(coefs(1,2)*v1(1)+coefs(2,2)*v1(2)+coefs(3,2)*v1(3));
dp2=(coefs(1,3)*v2(1)+coefs(2,3)*v2(2)+coefs(3,3)*v2(3));
dp3=(coefs(1,1)*v3(1)+coefs(2,1)*v3(2)+coefs(3,1)*v3(3));

res.dp=[dp1,dp2,dp3];
n=zeros(3,3);
n(:,1)=n1;
n(:,2)=n2;
n(:,3)=n3;
res.n=n;

function[res]=ev_plane_eq(x,eq)
res=eq(1)*x(1)+eq(2)*x(2)+eq(3)*x(3)+eq(4);

function[res]=is_inside_tetrahedron(x,tetr)
dscr=-1e-10;
if(tetr.deg~=0) res=0; return; end;
res=1;
for i=1:4
    if(ev_plane_eq(x,tetr.faces(:,i))<dscr)
        res=0; return;
    end;
end;

function[xr]=atb(x,sz,nDim)
xr=x;
for i=1:nDim
    if(x(i)<1) xr(i)=1; end;
    if(x(i)>sz(i)) xr(i)=sz(i);end;
end;
