function [m,m11,m12,m21,m22]=quadrant_mean(map,mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the mean of the whole map (except masked pix) and each quadrant
%Input:
%   -map
%   -mask
%Output:
%   -m:mean of whole map
%   -m11:mean of map(1:npix/2,1:npix/2)
%   -m12:mean of map(1:npix/2,npix/2+1:npix)
%   -m21:mean of map(npix/2+1:npix,1:npix/2)
%   -m22:mean of map(npix/2+1:npix,npix/2+1:npix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npix=size(map);npix=npix(1);
mm=map.*mask;m=mean(mm(find(mm)));
m11=mm(1:npix/2,1:npix/2);m11=mean(m11(find(m11)));
m12=mm(1:npix/2,npix/2+1:npix);m12=mean(m12(find(m12)));
m21=mm(npix/2+1:npix,1:npix/2);m21=mean(m21(find(m21)));
m22=mm(npix/2+1:npix,npix/2+1:npix);m22=mean(m22(find(m22)));
end
