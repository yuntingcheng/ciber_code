function [mask,clipmax,clipmin]=sigclip_mask(rawmap,rawmask,sig,iter_clip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function take a rawmap(unmasked or not well-masked),
%and raw mask (if no raw mask, just imput ones(n)). 
%and do the following two things:
%   -1. set nan to zero
%   -2. mask bad pixel for given iter_clip
%       and produce a new mask which mask out pixels in 
%       rawmask as well as crazy pixels
%
%Input:
%   -rawmap: orginal map (could also masked)
%   -rawmask: original mask
%   -sig: sigma to clip
%   -iter_clip: # of iteration to throw out crazy pixels
%Output:
%   -mask: mask inclue rawmask and crazy pixel
%   -clipmax:upper threashould for clipping
%   -clipmin:lower threashould for clipping
%   (clipmax, and clipmin are in the unit of map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set NaN to zero
rawmap(find(rawmap~=rawmap))=0;
rawmap(find(rawmap==inf))=0;
rawmap(find(rawmap==-inf))=0;
rawmask(find(rawmask~=rawmask))=0;
%% clip crazy pix and subtract mean
map1=rawmap.*rawmask;
[mask1,~,~,clipmax,clipmin]=sigma_clip(map1,sig,iter_clip);
mask=mask1.*rawmask;

return