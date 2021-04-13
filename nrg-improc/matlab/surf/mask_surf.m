% 06/25/2013: major updates
% - Input parameters are passed as structure from mask_surf_auto
% - Object mask can be specified as parameter
% - Computed object mask is saved.

% mask surface of the volume.
% params structure fields:
% root: root input file name
% opt: output more or less intermediate files. 0 - fewer (default, faster), 3 - most
% (slowest)
% method: 'blur','normfliter','coating'
% thresh: integer threshold between object and background for generating
% isosurface.
% grain: =1, default. <1, decrease % grid step. >1,
% increase % grid step.
% outroot: outfile name
% V0: input volume in original size and orientation.
% pixdim: vector with pixel sizes along x,y,z
% Vgfc: input volume used to calculate the binary mask (may be
% the same as V)
% img_exm: volume containing the exclusion mask (voxels that do
% not contribute to the volume mask and are prohibited to modify)
% img_head_mask: volume containing the head mask.
% inv: direction of casting rays to generate
% surface is from low to high coordinate (default is from high to low
% coordinate of vertical axis).
% xmin, xmax: outer and inner diagonal points of the bounding cuboid
% ver: can be 'x', 'y', or 'z'. denotes the anterior-posterior direction in
% the case of face

%function[R] = mask_surf(Vorig,Vgfc,Vexm, pixdim, root, vertical, bg_thresh, inverse, grid_step, method, opt)
%function[R] = mask_surf(Vorig,Vgfc,Vexm, pixdim, root, vertical, bg_thresh, inverse, grid_step, method, opt)
function[RR,msg] = mask_surf(params)

msg='unknown error';
nVols=uint16(size(params.V0,4));
xmin=params.xmin;
xmax=params.xmax;
pixdim=params.pixdim;

root=params.root;
vertical=params.ver;
bg_thresh=params.thresh;
inverse=params.inv;
grid_step=params.grain;
method=params.method;
opt=params.opt;

disp('computing object mask');
if (isempty(params.img_head_mask))
%computing mask
    if (bg_thresh==-1)
        %vertical partition into three volumes.
        [Vreor newdim] = reorient(params.Vgfc,vertical,1,inverse,pixdim,[0 0 0]);
        Mask=get_mask(Vreor);
        Mask=reorient(Mask,vertical,-1,inverse,newdim,[0 0 0]);
    else
        Mask=(params.Vgfc>=bg_thresh);
        disp(['Threshold: ' num2str(bg_thresh)]);
    end;
else
    Mask=params.img_head_mask;    
end;

if ~exist('newdim','var')
       newdim=pixdim;
end;

save_vol(Mask,params.V0_avw,[root '_obj_mask']);

Vorig=params.V0(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3),1:nVols);
Vexm=params.img_exm(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3));
Vgfc=params.Vgfc(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3));
Mask=Mask(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3));

[Mask0 ~]=reorient(Mask,vertical,1,inverse,pixdim,[0 0 0]);
[thickness, step] = calc_step(newdim, grid_step, size(Mask0));
pad=[step thickness];

[V newdim]=reorient(Vorig,vertical,1,inverse,pixdim,pad);
[Vgfc ~]=reorient(Vgfc,vertical,1,inverse,pixdim,pad);
[Vexm ~]=reorient(Vexm,vertical,1,inverse,pixdim,pad);
[Mask ~]=reorient(Mask,vertical,1,inverse,pixdim,pad);

Vexm=double(Vexm>0.5);
Mask=double(Mask>0.5);
%these are the default settings for 1 mm isotropic voxel.
%[thickness, step] = calc_step(newdim, grid_step, size(Mask));
disp(['thickness: ' num2str(thickness) ', step: ' num2str(step)]);
nm=12.0/norm(pixdim);
nm=max(1,ceil(nm));
disp(['nm: ' num2str(nm)]);
%morphologically close.
se=strel('disk',nm);
rx=(1/newdim(1));ry=(1/newdim(2));rz=(1/newdim(3));
%se=get_neib(max(1,round(2*rx)),max(1,round(2*ry)),max(1,round(2*rz)));

Mask=imopen(Mask,se);
Mask=imclose(Mask,se);

%ex_se=strel('disk',nm*5);
ex_se=get_neib(max(2,round(4*rx)),max(2,round(4*ry)),max(2,round(4*rz)));

Vexm=imdilate(Vexm,ex_se);
%Mask=Mask-Vexm;

%if(opt>0)
    saveVol(Mask,[root '_mask']);
%end;
if(opt>1)
    saveVol(V,[root '_orig']);
end;
%sz=size(V);
disp('computing boundary layer');
%calculate the volume mesh.
[faces, vertices, faces_b, vertices_b, faces_t,vertices_t,...
    lower_surf,upper_surf,orig_surf,mask_ex]=RectangularMesh(root,Mask,thickness,step,0.5,1.0,nm,opt);

%display generated mesh.

%f=figure('Visible','off');
%,'OuterPosition',[0,0,1024,768]);

disp ('saving the boundary layer');

%show fiture only if display is available.
if ~isdeployed
    ss=get(0,'ScreenSize');
    if(ss(3)>=100 || ss(4)>=100)    
        test_rect_mesh(faces, vertices, faces_b, vertices_b, faces_t,vertices_t);
    end
end

%debug
%return;

%print('-dpng',[ root '_mesh.png' ]);
%close(f);
%pause
disp ('masked surface rendering');

tic
Rb=[];
Rc=[];
Rnf=[];
Vmask=[];

nVols=size(Vorig,4);
for nv=1:nVols
    Frame=apply_exclude_mask(V(:,:,:,nv),mask_ex);
    %saveVol(Frame,'vol_exm');
    %saveVol(mask_ex,'mask_ex');
    %debug
    %continue;
    
    if(strcmp(method,'coating')~=0)
        disp ([ 'volume masking, frame ' num2str(nv)] );      
        Vmask=maskVol(Frame,lower_surf,orig_surf,upper_surf,step);
        if(nv==1 && opt>1)
            saveVol(Vmask,[root '_coatingMask']);
        end;
        disp (['finishing volume masking, frame ' num2str(nv)] );
        stats=get_seg_stats(Frame,Vmask,[1 2]);
        m=stats(1,1);
        Vrev1=m*Vmask+(1-Vmask).*Frame;
        Vrev=uint16(Vrev1.*(1-Vexm)+Vexm.*Frame);
        if(nv==1 && opt>0)
            saveVol(Vrev,[root '_coating']);
        end;
        %return to the original orientation.
        Vrev=reorient(Vrev,vertical,-1,inverse,newdim,pad);
        if(nv==1)
            R=zeros([size(Vrev) nVols]);
        end;
        R(:,:,:,nv)=Vrev;
        msg='ok';
    elseif(strcmp(method, 'blur')~=0)
        disp(['boundary layer blurring, frame ' num2str(nv)]);
        Vmask=maskVol(Frame,lower_surf,orig_surf,upper_surf,step);
        disp(['finalizing blur, frame ' num2str(nv)]);
        if(nv==1 && opt>1)
            saveVol(Vmask,[root '_blurMask']);
        end;
        Vbl=blur3(Frame,2*thickness);
        Vrev1=Vmask.*Vbl+(1-Vmask).*Frame;
        Vrev=uint16(Vrev1.*(1-Vexm)+Vexm.*Frame);
        if(nv==1 && opt>0)
            saveVol(Vrev,[root '_blur']);
        end;
        %return to the original orientation.
        Vrev=reorient(Vrev,vertical,-1,inverse,newdim,pad);
        if(nv==1)
            R=zeros([size(Vrev) nVols]);
        end;
        R(:,:,:,nv)=Vrev;
        msg='ok';
    elseif(strcmp(method, 'normfilter')~=0)
        %project forward.
        disp(['projecting face, frame ' num2str(nv)]);
        Vres=projectVol(Frame,Frame,orig_surf,lower_surf,upper_surf,thickness,step,'direct',opt);
        if(opt>0)
            saveVol(Vres,[root '_face']);
        end;
        %mask volume.
        disp('lowpass filtering');
        V1=blur3_thin(Vres,thickness);
        if(opt>0)
            saveVol(uint16(V1),[root '_face_blur']);
        end;
        %project backward.
        disp('back projection');
        [Vrev1,Vmask]=projectVol(V1,Frame,orig_surf,lower_surf,upper_surf,thickness,step,'reverse',opt);
        Vrev=uint16(Vrev1.*(1-Vexm)+Vexm.*Frame);
        disp(['finalizing normalized filtering, frame ' num2str(nv)]);
        %return to the original orientation.
        if(opt>0)
            saveVol(uint16(Vrev),[root '_normfilter']);
        end;
        Vrev=reorient(Vrev,vertical,-1,inverse,newdim,pad);
    %   saveVol(uint16(Vrev),[root '_normfilter']);    
        if(nv==1)
            R=zeros([size(Vrev) nVols]);
        end;
        R(:,:,:,nv)=Vrev;
        msg='ok';
    end;    
    if(isempty(Vmask)) display ('invalid method requested, exiting'); exit; end;
    if (nv==1 && opt>0)
        saveVol(Vmask,[root '_layer']);
    end;    
end;    
%Vmask=reorient(Vmask,vertical,-1,inverse);
RR=params.V0;
RR(xmin(1):xmax(1),xmin(2):xmax(2),xmin(3):xmax(3),:)=R;
disp('surface masking done');
toc;
return;



function[Mask]=get_mask(V)
Mask=zeros(size(V));
sz=size(V);
h=sz(2);
h0=uint16(h/3); h1=uint16(2*h0); h2=h;
t1=select_threshold(V(:,1:h0,:));
t2=select_threshold(V(:,h0:h1,:));
t3=select_threshold(V(:,h1:h2,:));
Mask(:,1:h0,:)=V(:,1:h0,:)>t1;
Mask(:,h0:h1,:)=V(:,h0:h1,:)>t2;
Mask(:,h1:h2,:)=V(:,h1:h2,:)>t3;
disp (['thresholds: ' num2str(t1) ' ,' num2str(t2) ', ' num2str(t3)]);

function[thickness,step]=calc_step(pixdim,grid_step,dims)
%factor=grid_step;
%pd=pixdim;
%if(grid_step<0.8) 
%    pd=ones(3)*max(pd(:));
%elseif(grid_step>1.2)
%    pd=ones(3,1)*min(pd(:));
%end
%if(pd(1)<=0 || pd(2)<=0 || pd(3)<=0) return; end;
%diag=sqrt(sum(pd.*pd));
%step=double(max(4,round((factor*10*sqrt(3))/diag)));
zlim=[1 floor((dims(3)-1)/2)];
xlim=[1 floor((dims(1)-1)/2)];
ylim=[1 floor((dims(2)-1)/2)];

thickness=min(zlim(2),max(zlim(1),round(grid_step*8.0/pixdim(3))));

wx=round(grid_step*15/pixdim(1));
wy=round(grid_step*15/pixdim(2));
step(1)=min(xlim(2),max(xlim(1),wx));
step(2)=min(ylim(2),max(ylim(1),wy));

%thickness=double(max(2,round((factor*6*sqrt(3))/diag)));

function[res newdim]=reorient(V, vertical, dir, inverse, pixdim, pad)

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
       % res=reflect(res,3);
        res=flipdim(res,3);
    end;
    %pad with zeroes
    res=padarray(res,pad,'replicate');
else %back
    d=size(V);
    %unpad
    res=V(pad(1)+1:d(1)-pad(1),pad(2)+1:d(2)-pad(2),pad(3)+1:d(3)-pad(3),:);
    res=permute(res,transp);
    if(inverse==1)
        if(strcmp(vertical,'x'))
           % res=reflect(res,1);
            res=flipdim(res,1);
        elseif(strcmp(vertical,'y'))
           % res=reflect(res,2);
            res=flipdim(res,2);
        else 
           % res=reflect(res,3);
            res=flipdim(res,3);
        end
    end
end;
function[M] = get_neib(rx,ry,rz)
M=zeros(2*rx,2*ry,2*rz);
cx=rx+.5;cy=ry+.5;cz=rz+.5;
rx2=1/(rx*rx);ry2=1/(ry*ry);rz2=1/(rz*rz);
for z=1:2*rz+1
    for y=1:2*ry+1
        for x=1:2*rx+1
            if (rx2*(x-cx)^2+ry2*(y-cy)^2+rz2*(z-cz)^2 <= 1)
                M(x,y,z)=1;
            end;
        end;
    end;
end

function[V1]=apply_exclude_mask(V,mask_ex)
    sz=size(V);
    %compute median    
    med=prctile(V(mask_ex>0),50);
    V1=V; V1(mask_ex>0)=med;
