function [sm,fillmap]=fillpadsmooth(map,mask,sig,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given a map and its mask, apply Gaussian smooth with sigma
%after padding it to an area twice the size to avoid edge effects.
%
%Input:
%  -map
%  -mask
%  -sig: sigma of the Gaussian filter
%  -hsize_r: hsize of Gaussian filter w.r.t sig, default:5
%   (set hsize=round(sig*hsize_r))
%
%Output:
%  -sm: smoothed map 
%  -fillmap: fillmap return by unholy_map(map,mask,100)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('map',@isnumeric);
  p.addRequired('mask',@isnumeric);
  p.addRequired('sig',@isnumeric);
  p.addOptional('hsize_r',5,@isnumeric);
  
  p.parse(map,mask,sig,varargin{:});

  map     = p.Results.map;
  mask     = p.Results.mask;
  sig = p.Results.sig;
  hsize_r    = p.Results.hsize_r;
  
  clear p varargin;


s = size(map);cx=s(1);dx = round(cx/2);

bigun = zeros(s(1)*2);

fillmap = unholy_map(map,mask,100);
bigun(cx-dx:cx+dx-1,cx-dx:cx+dx-1)=fillmap;

gauss = fspecial('Gaussian',round(sig*hsize_r),sig);
zers = find(bigun == 0);

bigun(zers) = mean(fillmap(:));
bigsm = imfilter(bigun,gauss);
sm = bigsm(cx-dx:cx+dx-1,cx-dx:cx+dx-1);

end