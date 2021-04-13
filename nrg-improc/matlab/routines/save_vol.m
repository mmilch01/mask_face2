% Author: 	Mikhail Milchenko, mmilch@npg.wustl.edu, Washington University School of Medicine. 
% Warranty:	No express or implied warranty, or fitness of this code for specific purpose is assumed. 
% License:	You are free to use and re-distribute this code in any form, preserving the above notice.

function []=save_vol(V,avw,fname)
sz=size(V);
vavw=avw;
vavw.hdr.dime.dim(1)=4;
vavw.hdr.dime.dim(2)=sz(1);
vavw.hdr.dime.dim(3)=sz(2);
vavw.hdr.dime.dim(4)=sz(3);
if(ndims(V)>3)
    vavw.hdr.dime.dim(5)=sz(4);
else
    vavw.hdr.dime.dim(5)=1;
end;
%vavw.img=V;
vavw.hdr.dime.datatype=16;
vavw.hdr.dime.bitpix = int16(32);
vavw.img=single(V);
avw_img_write(vavw,fname,0);