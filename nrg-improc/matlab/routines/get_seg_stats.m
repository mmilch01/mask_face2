% Returns sample statistics (2xnClasses) of tissue classes
% I - source volume, Iseg-label volume, Seg-label array

% Author: 	Mikhail Milchenko, mmilch@npg.wustl.edu, Washington University School of Medicine. 
% Warranty:	No express or implied warranty, or fitness of this code for specific purpose is assumed. 
% License:	You are free to use and re-distribute this code in any form, preserving the above notice.

function[out]= get_seg_stats(I,Iseg,Seg)
nClasses=max(size(Seg));
out=zeros(2,nClasses);
for i=Seg(1):Seg(nClasses)
    t=min(max(i-1,Iseg),i+1)-i;
    t=-t.*t+1;
    n=sum(sum(sum(t)));
    t1=t.*I;
    s1=sum(sum(sum(t1)))/n;
    s2=t.*(t1-s1);
    out(1,i)=s1;
    out(2,i)=sum(sum(sum(s2.*s2)))/n;
end;
return;