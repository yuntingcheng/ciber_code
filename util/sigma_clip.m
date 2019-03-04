function [mask,frac1,frac2,clipmax,clipmin]=sigma_clip(maskmap,sig,iter_clip)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input a maskedmap, and mask out bad pixels in the maskmap.
%Output:
%   - mask: new mask to mask out crazy pixels
%     Note!!only crazy pixels are masked, the
%     original mask does not included!!!
%   -frac1:fraction of non-zero/non-nan pixel in input maskmap
%   -frac2:fraction of pixels not masked of output mask
%   -clipmax:upper threashould for clipping
%   -clipmin:lower threashould for clipping
%   (clipmax, and clipmin are in the unit of map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=maskmap(:);
b=b(~~b); %only take non-zero values

frac1=numel(b)/numel(maskmap);
  for i=1:iter_clip
      b = b(b < nanmedian(b(:))+sig*std(b(:)) ... 
          & b > nanmedian(b(:))-sig*std(b(:)));
  end
  
%  b = b(b < nanmedian(b(:))+5*std(b(:)) ...
%     & b > nanmedian(b(:))-5*std(b(:)));
  
mask=maskmap<max(b) & maskmap>min(b);
mask=double(mask);
clipmax=max(b);clipmin=min(b);
b1=b(~~b);
frac2=numel(b1)/1024.^2;
return
