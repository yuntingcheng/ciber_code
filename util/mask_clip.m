function [mask,frac,clipmax,clipmin]=mask_clip(map,iter_clip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function produce a mask to mask out
%the crazy pixel. The masked part is 0,
%and the not masking part is 1. Multiply 
%the map by this mask can get the masked map.
%
%Output:
%   - iter_clip: # of iteration used for clip
%   (same method as imageclip.m)
%   -frac:fraction of pixels non-masked
%   -clipmax:upper threashould for clipping
%   -clipmin:lower threashould for clipping
%   (clipmax, and clipmin are in the unit of map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b=map(:);
  for i=1:iter_clip
      b = b(b < nanmedian(b(:))+5*std(b(:)) ... 
          & b > nanmedian(b(:))-5*std(b(:)));
  end
b = b(b < nanmedian(b(:))+5*std(b(:)) ...
      & b > nanmedian(b(:))-5*std(b(:)));
clipmax=max(b);clipmin=min(b);
  
frac=numel(b)/numel(map);
mask=map<max(b) & map>min(b);

return
