%06/25/2013: major update of interfaces
% - Variable lengh input
% - all input parameters are computed in a separate function and passed to
% mask_surf
% - Can take manual ROI's
% - Saves computed ROI's; can input previously saved ROI's
% - Saves computed head mask; can input previously computed head mask
% - saves computed directional normal to surface; can input directional
% normal
% - file in ROI format is saved under specified label
% 02/29/2013: Van Gogh update.
% 06/18/2013: Reformat the input.
% 12/27/2013: added accounting for highly anisotropic images; ability to
% save snapshots only.

function[]=mask_surf_auto(root,varargin)
if isdeployed
    v1=convert_num_args(varargin);
else
    v1=varargin;
end

p=inputParser;
%input image root name.
p.addRequired('root',@ischar);
%ROI coordinate file name.
p.addParamValue('coords',[],@ischar);
%output image root name (default is <root>_full_<method>)
p.addParamValue('outroot',[],@ischar); 
%exclusion mask image file name
p.addParamValue('exmask',[],@ischar);
%masking method, normfilter is default.
p.addParamValue('method','normfilter',@(x)strcmpi(x,'blur') | strcmpi(x,'normfilter') | strcmpi(x,'coating'));
%object-background intensity threshold, -1 for auto threshold.
p.addParamValue('thresh',-1,@(x)isnumeric(x));
%mesh step factor, 1 is default, smaller/bigger will result in less/more
%defacing
p.addParamValue('grain',1,@(x)isnumeric(x) && (x>=.3 && x<=3));
%output more or less intermediate files. 0 - fewer (default, faster), 3 - most
%(slowest)
p.addParamValue('optimization',0,@(x)isnumeric(x) && (x==0 || x==1 || x==2 || x==3));
%pre-calculated head mask image can be supplied. 
p.addParamValue('head_mask',[],@ischar);
% all of pmin,pmax,vertical, must be either defined or undefined
% simultaneously.
%outer 3D ROI pinnacle, furthest from the face surface
p.addParamValue('p_out',[],@(x)isnumeric(x));
%inner 3D ROI pinnacle, closest to the face surface.
p.addParamValue('p_in',[],@(x)isnumeric(x)); 
%axis that is most normal to the face surface.
p.addParamValue('vertical',[],@(x)isnumeric(x) && (x==0 || x==1 || x==2));
%instead of vertical, specify a directed normal (e.g. [0 -1 0] would mean
%vertical=1 and inverse=1
p.addParamValue('normal',[],@(x)isnumeric(x)); 
%ROI label file name.
p.addParamValue('roi_label',[],@ischar);
%Option to only output 3D snapshot.
p.addParamValue('snapshotOnly',[],@(x)isnumeric(x));

p.parse(root,v1{:});
disp('Input arguments:'); p.Results
%initialize parameters.
[params msg]=init_params(p);
if isempty(params) disp(msg); return; end
disp('Calculated parameters:'); params

if (params.snapshotOnly==0)
    %save ROI file. 
    save_ROI(params, strcat(params.roi_label,'.roi'));

 %   snapshot(params.V0,params.pixdim,params.ver,params.inv,[params.root '_orig']);
 %   snapshot3D(params.V0,params.pixdim,params.ver,params.inv,[params.root '_orig_surf']);

    %run surface masking using the pre-calculated parameters.
    [R msg] = mask_surf(params);
    if (isempty(R)) disp(msg); return; end
    avw=avw_img_read(strcat(params.root,'.img'),0);
    avw.img=R;
    %save the result.
    save_vol(uint16(R),avw,params.outroot);
    %save QC snapshots.
    snapshot(R(:,:,:,1),params.pixdim,params.ver,params.inv,[params.outroot '_' params.method '.png']);
    snapshot3D(R(:,:,:,1),params.pixdim,params.ver,params.inv,[params.outroot '_' params.method '_surf']);
else
    R=params.V0;
    %save snapshots.
    snapshot(R(:,:,:,1),params.pixdim,params.ver,params.inv,[params.outroot '_' params.method '.png']);
    snapshot3D(R(:,:,:,1),params.pixdim,params.ver,params.inv,[params.outroot '_' params.method '_surf']);
end;
return;

%initialize all parameters needed to run the mask_surf.
function [params,message] = init_params(p)
message='ok';

%first, init simple params.
params.root=p.Results.root;
params.opt=p.Results.optimization;
params.method=p.Results.method;
params.thresh=p.Results.thresh;
params.grain=p.Results.grain;

if isempty(p.Results.outroot)
    params.outroot=strcat(params.root,'_full_',params.method);
else
    params.outroot=p.Results.outroot;
end;

if isempty(p.Results.roi_label)
    params.roi_label=params.root;
elseif isempty(p.Results.outroot)
    params.outroot=p.Results.roi_label;
end;

if isempty(p.Results.snapshotOnly)
    params.snapshotOnly=0;
else
    params.snapshotOnly=1;
end

params.roi_label=p.Results.roi_label;

%load images.
im=strcat(params.root,'.img');

avw=avw_img_read(im,0);
params.V0=avw.img;
params.V0_avw=avw;
params.pixdim=double(avw.hdr.dime.pixdim(2:4));
disp(['pixdim: ' num2str(params.pixdim)]);
params.Vgfc=params.V0;


if exist (strcat(p.Results.exmask,'.img'),'file')
    avwem=avw_img_read(strcat(p.Results.exmask,'.img'),0);
    params.img_exm=fixbet(avwem.img);
else 
    params.img_exm=params.V0*0;
end      
if exist (strcat(p.Results.head_mask,'.img'), 'file')
    avwhm=avw_img_read(strcat(p.Results.head_mask,'.img'),0);
    params.img_head_mask=avwhm.img;
else 
    params.img_head_mask=[];
end

%calculate the ROI boundaries.
pts=p.Results.coords;
ver=p.Results.vertical;
xmin=p.Results.p_out;
xmax=p.Results.p_in;
normal=p.Results.normal;
%manually specified ROI
if ~isempty(xmin) && ~isempty(xmax) && (~isempty(ver) || ~isempty(normal))
    sz=size(params.V0);
    xmin=uint16(boundcheck(p.Results.p_out,params.V0));
    xmax=uint16(boundcheck(p.Results.p_in,params.V0));
%    if isempty(normal) %flip y coordinate from the manual ROI.
%        xmin=[xmin(1) sz(2)-xmin(2) xmin(3)];
%        xmax=[xmax(1) sz(2)-xmax(2) xmax(3)];
%    end
    d=xmax(ver+1)-xmin(ver+1);
    if(d>0) inv=1; else inv=0; end;
    if(~isempty(normal)) %
        [tmp ind]=max(abs(normal));
        ver=ind-1;
        inv=abs(min(0,sign(normal(ver+1))));
    end
    params.inv=inv;    
    
    switch ver
        case 0
            params.ver='x';
        case 1
            params.ver='y';
        case 2 
            params.ver='z';
    end
    params.xmin=min(xmin,xmax);
    params.xmax=max(xmin,xmax);
%read the pre-calculated ROI file.    
elseif exist (pts,'file') 
    coord=dlmread(pts);
    dims=size(params.V0(:,:,:,1));
    [x0,xmin,xmax,d]=ROI_coord(params.V0(:,:,:,1),coord);
    b=(xmin==xmax);
    c=(isnan(x0)+isnan(xmin)+isnan(xmax));
    if(sum(b)>0 || sum(c)>0)
        params=[];
        message='RESULTING ROI IS TOO SMALL FOR DEFACING';
        return;
    end;
    if(d(1)>0) inv=0; else inv=1; end;
    switch d(2)
        case 1
            ver='x';
            if (inv==1) xmin(1)=1; else xmax(1)=dims(1); end;
        case 2
            ver='y';
            if (inv==1) xmin(2)=1; else xmax(2)=dims(2); end;
        case 3
            ver='z';
            if (inv==1) xmin(3)=1; else xmax(3)=dims(3); end;
        otherwise
            params=[];
            message='Cannot calculate surface ROI, exiting';
            return;
    end;
    params.inv=inv;
    params.ver=ver;
    %increase the ROI to compensate for partial segment effect.
    nm=round(6*params.pixdim);
    params.xmin=int16(max(xmin-nm,ones(1,3)));
    params.xmax=int16(min(xmax+nm,dims));
else
    message='Cannot calculate surface ROI, exiting';
    params=[];
    return;   
end
%save normal.
if isempty(normal)
    switch params.ver
        case 'x'
            normal=[1 0 0];
        case 'y'
            normal=[0 1 0];
        case 'z'
            normal=[0 0 1];
    end;
    if (params.inv==1)
        normal=normal*-1;
    end;
    params.normal=normal;
end

function []=save_ROI(params, fname)
switch params.ver
    case 'x'
        n=[1 0 0];
    case 'y'
        n=[0 1 0];
    case 'z'
        n=[0 0 1];
end
if (params.inv==1) n=n*-1; end;
out=double([double(params.xmax); double(params.xmin); n]);
save(fname,'out','-ASCII');

    
function[x0,xmin,xmax,d]=ROI_coord(I0,coord)
orig = coord(1,:);
x=coord(2:21,:);
y=coord(22:41,:);
z=coord(42:61,:);
[~, x00]=getROI2D(I0,orig);
[x10 x11]=getROI2D(I0,x);
[x20 x21]=getROI2D(I0,y);
[x30 x31]=getROI2D(I0,z);
%v1=x1-x0; %v1=v1/norm(v1,'fro');
v20=x11-x00;
v21=x10-x00;
v2=v21;
if(norm(v20,'fro')>norm(v21,'fro')) v2=v20; end;
v2=v2/norm(v2,'fro');
%v3=x3-x0; %v3=v3/norm(v3,'fro');
e1=[1 0 0];
e2=[0 1 0];
e3=[0 0 1];
d=[dot(v2,e1) dot(v2,e2) dot(v2,e3)];
x0=x00;
xmin=min([x10; x11; x20; x21; x30; x31]);
xmax=max([x10; x11; x20; x21; x30; x31]);
[~, j]=max(abs(d));
val=d(j);
d=[val j];

function[] = snapshot(V,pixdim, ver, inv, name)
%f=figure('Visible','off','OuterPosition',[0,0,1024,768]);
[Vr newdim]=reorient(V(:,:,:,1),ver,1,inv,pixdim);
f=figure('Visible','on','OuterPosition',[0,0,1280,900]);
Vr=permute(Vr,[3 1 2]);
Vr=flipdim(Vr,1);
Vr=flipdim(Vr,3);
dispvol(Vr,[],12);
%print('-dpng', name);
sv(f,name);
close(f);

function[] = snapshot3D(V,pixdim,ver,inv,name)
%memory;
[Vr newdim]=reorient(V,ver,1,inv,pixdim);
thresh=select_threshold(V);
%f=figure('Visible','off','OuterPosition',[0,0,1024,768]);
%f=figure;
I=dispvol3D(Vr(:,:,:,1),newdim,thresh);
imwrite(I,[name '.png'],'png');
%sv(f,name);
%print('-dpng',name);
%close(f);

function sv(f, name)
%im=frame2im(getframe(f));
%imwrite(im,[name '.png'],'png');
%export_fig(name);
%opengl software;
saveas(f,name,'png');
%set(gcf,'Renderer','zbuffer');
%print('-dpng',name);
%hardcopy(gcf, '-Dopengl', '-r300');

function[mn, mx]=getROI2D(V,arr)
if(size(arr,1)==1) 
    mn=arr; mx=arr; return;
end;    
sz=size(V);
mn=max([1 1 1], min(arr));
mx=min(sz,max(arr));

function[res newdim]=reorient(V, vertical, dir, inverse, pixdim)
is4d=(size(V,4)>1);
if(strcmp(vertical,'x'))
   transp=[3 2 1];
elseif(strcmp(vertical,'y'))
   transp=[1 3 2];
else 
   transp=[1 2 3];
end;
newdim=[pixdim(transp(1)) pixdim(transp(2)) pixdim(transp(3))];
if (is4d)
    transp=[transp 4];
end;

if(dir==1) %forward
    res=permute(V,transp);
    if(inverse==1)
        res=flipdim(res,3);
       % res=reflect(res,3);
    end;
else %back
    res=permute(V,transp);
    if(inverse==1)
        if(strcmp(vertical,'x'))
%            res=reflect(res,1);
             res=flipdim(res,1);
        elseif(strcmp(vertical,'y'))
%            res=reflect(res,2);
            res=flipdim(res,2);
        else 
%            res=reflect(res,3);
            res=flipdim(res,3);
        end
    end
end;

function [y]=boundcheck(x,V)
dim=size(V(:,:,:,1));
y=min(max([1,1,1],x),dim);

function[Vm]=fixbet(V)
dx=size(V,1);
dy=size(V,2);
dz=size(V,3);
x1=V(1,:,:);
x2=V(dx,:,:);
y1=V(:,1,:);
y2=V(:,dy,:);
z1=V(:,:,1);
z2=V(:,:,dz);
a=2*(dx*dy+dx*dz+dy*dz);
a1=sum(x1(:))+sum(x2(:))+sum(y1(:))+sum(y2(:))+sum(z1(:))+sum(z2(:));
if( a1/a > .5)
    Vm=ones(size(V))-V;
else
    Vm=V;
end;

function [out]=convert_num_args(v)
numargs={'thresh','grain','optimization','p_out','p_in','vertical','normal','snapshotOnly'};
sz=size(v,2)
for i=1:sz
    for j=1:size(numargs,2)
        if strcmp(v(i),numargs(j))==1
            m=str2num(char(v(i+1))); sm=size(m);
            v(i+1)=mat2cell(m,sm(1),sm(2));
            i=i+1;
        end
    end    
end
out=v;
