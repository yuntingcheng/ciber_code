function [map,mask,meanmap,frac1,frac2,clipmax,clipmin]=...
        get_skymap(rawmap,rawmask,sig,iter_clip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function take a rawmap(unmasked or not well-masked),
%and raw mask (if no raw mask, just imput ones(n)). 
%and do the following two things:
%   -1. set nan to zero
%   -2. throw away bad pixel for given iter_clip
%       and produce a new mask which mask out pixels in 
%       rawmask as well as crazy pixels
%   -3. subtract mean to make the final map zero-meaned
%
%Input:
%   -rawmap: orginal map (could also masked)
%   -rawmask: original mask
%   -iter_clip: # of iteration to throw out crazy pixels
%Output:
%   -map: masked-map which mask out rawmask and crazy pixel
%         and mean-substracted
%   -mask: mask inclue rawmask and crazy pixel
%   -frac1:fraction of non-zero/nan pixel in input maskmap
%   -frac2:fraction of pixels not masked of output mask
%   -clipmax:upper threashould for clipping
%   -clipmin:lower threashould for clipping
%   (clipmax, and clipmin are in the unit of map)
%   Note: clipmax and clipmin are also mean subtracted!!
%         so it is not same as mask_clip1 output!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set NaN to zero
rawmap(find(rawmap~=rawmap))=0;
rawmap(find(rawmap==inf))=0;
rawmap(find(rawmap==-inf))=0;
rawmask(find(rawmask~=rawmask))=0;
%% clip crazy pix and subtract mean
map1=rawmap.*rawmask;
[mask1,frac1,frac2,clipmax,clipmin]=sigma_clip(map1,sig,iter_clip);
map2=map1.*mask1;
meanmap=mean(map2(find(map2)));
mask=mask1.*rawmask;
map=(map2-meanmap).*mask;

clipmax=clipmax-meanmap;
clipmin=clipmin-meanmap;
return