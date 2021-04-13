function[] = dispvol(V0, M, nIm)
%V=permute(V0,[2 1 3]);
V=V0(:,:,:,1);
colormap(gray(256));

if ~exist('M','var')
    M=[]
end
if isempty(M)
    if size(unique(V(:)),1)>99
        M=prctile(V(:),98);
    else
        M=max(V(:));
    end
else
    colormap(hot(256));
    for i=1:numel(V)
        if(V(i)>M) V(i)=M; end;
    end;
end;

if ~exist ('nIm','var')
    nIm=[];
end;
if isempty(nIm)
    nIm=20;
end;

m=min(min(min(V)));

[dx dy dz] = size(V);

nRows=4;
nCols=6;
if(isempty(nIm))
    nIm=min(dz,24);
end;
if(nIm<=20)
    if(nIm>15)
        nCols=5; nRows=4;
    elseif (nIm>12)
        nCols=5; nRows=3;
    elseif (nIm>9)
        nCols=4; nRows=3;
    elseif (nIm>8)
        nCols=3; nRows=3;
    elseif (nIm>6)
        nCols=4; nRows=2;
    elseif (nIm>4)
        nCols=3; nRows=2;
    elseif (nIm>=3)
        nCols=2; nRows=2;
    elseif (nIm==2)
        nCols=2; nRows=1;
    else
        nCols=1; nRows=1;
    end;
end;
step = dz/nIm;
if(dz<nIm) step=1; end;

slice=-1;
ind=0;

%colormap(hot(256));
for i=1:nIm
    prev_slice=slice;
    if (i==1) slice=1;
    elseif (i==nIm) slice=dz;
    else slice=int16(i*step);
    end;
    if(prev_slice==slice) continue; end;
    ind=ind+1;    
    subplot(nRows,nCols,ind);
    imagesc(V(:,:,slice),[m M]),title(['slice ' num2str(slice)]);
%    imshow(V(:,:,slice),[m M]),title(['slice ' num2str(slice)]);
end;
colorbar;

function num_greylevels(V)

