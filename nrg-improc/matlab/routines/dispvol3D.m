function[Z] = dispvol3D(V,pixdim,thresh,angles)
%this line is for the batch offscreen figure saving.
%opengl software;

%subplot(1,2,1);
%Ds=smooth3(resample3D(V,pixdim,1));
if ~exist ('angles','var')
   angles=[-45 0 0];
end
dsz=0.8;
Ds=smooth3(resample3D(V,pixdim,dsz));
iso=isosurface(Ds,thresh);

%render the volume.
%Z=render(iso.vertices,iso.faces,[180 -80 0]);
%Z=render(iso.vertices,iso.faces,[90 -45 -90]);
%Z=render(iso.vertices,iso.faces,[135 180 0]);
Z=abs(render(iso.vertices,iso.faces,angles));
%imshow(Z,[0,1]);
return;

daspect(pixdim);
hiso=patch(iso,'FaceColor',[1,.75,.65],'EdgeColor','none');
%set(gcf,'Renderer','OpenGL'); lighting gouraud
%isonormals(Ds,hiso);
set(gcf,'Renderer','zbuffer'); lighting phong
view(-80,10);
axis tight;
%set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);
lightangle(-140,20);


%uncomment the following to add sagittal view.
%subplot(1,2,2);
%hiso=patch(iso,'FaceColor',[1,.75,.65],'EdgeColor','none');
% set(gcf,'Renderer','zbuffer'); lighting phong
% isonormals(Ds,hiso);
% set(hiso,'SpecularColorReflectance',0,'SpecularExponent',50);
% view(-80,10);
% lightangle(-140,20);
% view(-20,10);
% lightangle(-50,30);
% daspect(pixdim);
% axis tight;
% 

function[X]=resample3D(V,pixdim,sz)
osz=size(V);
newdims=osz.*pixdim/sz;
maxdims=[192 192 192];
exceed=sum(newdims>maxdims);
if(exceed>0)
    exratio=max(newdims./maxdims);
    sz=sz*exratio;
end
x=pixdim(1):pixdim(1):osz(1)*pixdim(1);
y=pixdim(2):pixdim(2):osz(2)*pixdim(2);
z=pixdim(3):pixdim(3):osz(3)*pixdim(3);
[yi,xi,zi]=meshgrid(1:sz:osz(1)*pixdim(1),1:sz:osz(2)*pixdim(2),1:sz:osz(3)*pixdim(3));

%X=interp3(x,y,z,V,xi,yi,zi,'linear');

X=interp3(y,x,z,V,xi,yi,zi,'linear');
X=permute(X,[2 1 3]);
return;


function[Z]=render(verts,faces,ang)

%initialize the angle calculation.
ang=ang*pi()/180;
fx=ang(1);fy=ang(2);fz=ang(3);
Rx=[[1 0 0];[0 cos(fx) sin(fx)];[0 -sin(fx) cos(fx)]];
Ry=[[cos(fy) 0 sin(fy)];[0 1 0];[-sin(fy) 0 cos(fy)]];
Rz=[[cos(fz) sin(fz) 0];[-sin(fz) cos(fz) 0];[0 0 1]];
R=Rx*Ry*Rz;

nv=verts*R;
shift=zeros(size(verts));
shift(:,1)=-min(nv(:,1))+1;
shift(:,2)=-min(nv(:,2))+1;
verts=verts*R+shift;

%initialize resulting matrices.
ymax=ceil(max(verts(:,2)));
xmax=ceil(max(verts(:,1)));

Z=zeros(xmax,ymax);
Zi=zeros(xmax,ymax);
nFaces=size(faces,1);
for i=1:nFaces
    face=faces(i,:);
    x=verts(face,:);
    x1=x(1,:);
    d1=x(2,:)-x(1,:);
    d2=x(3,:)-x(1,:);
    n1=norm(d1);
    n2=norm(d2);
    if(n1>n2)
        x2=x(2,:);
        x3=x(3,:);
    else        
        n=n1;
        n1=n2;
        n2=n;
        x2=x(3,:);
        x3=x(2,:);
    end;
    %check if there are intersections with the z plane.
    if (Z(round(x1(1)),round(x1(2)))>x1(3) && Z(round(x2(1)),round(x2(2)))>x2(3) && Z(round(x3(1)),round(x3(2)))>x3(3))
        continue;
    end;
    d=ceil(n1);
    dt=1/(d+1);
    for t=0:dt:n1
        u1=x1+n1*t;
        u2=x2+n2*t;
        dd=u2-u1;
        nn=norm(dd(1:2));
        dtt=1/(ceil(nn)+1);
        for tt=0:dtt:nn
            u=u1+dd*tt;
            ru=round(u);
            if(ru(1)<1 || ru(2)<1 || ru(1)>xmax || ru(2)>ymax) continue; end;
%            if(sum(ru(1:2)<[1 1])>0 || sum(ru(1:2)>[xmax ymax])>0) continue; end;
            if(Z(ru(1),ru(2))<u(3))
                Z(ru(1),ru(2))=u(3);
                Zi(ru(1),ru(2))=i;
            end;
 %          break;
        end;
        break;
    end;
end;

%render all visible triangles.
%light vector.
light=[0 0 -1];
light=light/norm(light);
maxval=255;
for x=1:xmax
    for y=1:ymax
        i=Zi(x,y);
        if(i==0) continue; end;
        face=verts(faces(i,:),:);
        nrm=cross(face(2,:)-face(1,:),face(3,:)-face(1,:));
        n=norm(nrm);
        if(n==0) 
            Z(x,y)=1;
        else
            Z(x,y)=dot(nrm,light)/n;
        end;
    end;
end;