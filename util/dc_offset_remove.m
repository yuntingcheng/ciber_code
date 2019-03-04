function [map]=dc_offset_remove(map,mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function mitigate the DC offset in dark map.
%This function gives a dc constant to each quadrant such that
%they have the same mean as the original map.
%Input:
%   -map
%   -mask
%Output:
%   -cCl:Cl after mkk
%   -wCl:Cl before mkk
%   -Cl2d: 2D PS before mkk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
npix=size(map);npix=npix(1);
[m,m11,m12,m21,m22]=quadrant_mean(map,mask);
map(1:npix/2,1:npix/2)=map(1:npix/2,1:npix/2)+m-m11;
map(1:npix/2,npix/2+1:npix)=map(1:npix/2,npix/2+1:npix)+m-m12;
map(npix/2+1:npix,1:npix/2)=map(npix/2+1:npix,1:npix/2)+m-m21;
map(npix/2+1:npix,npix/2+1:npix)=map(npix/2+1:npix,npix/2+1:npix)+m-m22;
end
