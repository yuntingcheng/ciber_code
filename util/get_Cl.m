function [cCl,wCl,Cl2d,l,binl,dCl]=get_Cl(map,mask,Mkk,beam,weight,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given the map and mask and the mkk of the mask, get 
%1DPS (wCl) and 2DPS(Cl2d) before Mkk, and after Mkk (cCl).
%Mask NaN, inf, -inf pixels before taking PS.
%Input:
%   -map
%   -mask
%   -Mkk
%   -beam[arcsec]
%   -weight: Fourier weight in 2DPS
%   -clipCl2d(optional): 1 or 0, doing sigma clip in Cl2d,default:0
%Output:
%   -cCl:Cl after mkk
%   -wCl:Cl before mkk
%   -Cl2d: 2D PS before mkk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Parse data
  p = inputParser;
  
  p.addRequired('map',@isnumeric);
  p.addRequired('mask',@isnumeric);
  p.addRequired('Mkk',@isnumeric);
  p.addRequired('beam',@isnumeric);
  p.addRequired('weight',@isnumeric);
  p.addOptional('clipCl2d',0,@isnumeric);
  
  p.parse(map,mask,Mkk,beam,weight,varargin{:});

  map     = p.Results.map;
  mask     = p.Results.mask;
  Mkk     = p.Results.Mkk;
  beam     = p.Results.beam;
  weight     = p.Results.weight;
  clipCl2d     = p.Results.clipCl2d;
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
map(find(map~=map))=0;
map(find(map==inf))=0;
map(find(map==-inf))=0;

map=map.*mask;
map=map-mean(map(find(map)));map=map.*mask;
[~,l,~,~,binl,~,Cl2d] = get_angular_spec(map,map,beam);

if clipCl2d==0
[wCl]=Cl_from_Cl2d(Cl2d,beam,'w',weight);
else
[wCl]=Cl_from_Cl2d_clip(Cl2d,beam,'w',weight);
end  

tcl = find(wCl~=0.0 & ~isnan(wCl));
Mkkp = Mkk(tcl,tcl);
cCl = zeros(1,numel(tcl));
cCl(tcl) = (inv(Mkkp)*wCl(tcl)')';

[dCl]=dCl_Knox(cCl,binl,beam);

l=l(find(cCl));wCl=wCl(find(cCl));
dCl=dCl(find(cCl));cCl=cCl(find(cCl));
end
