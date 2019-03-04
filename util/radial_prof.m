function profile = radial_prof(map,mask,cenx,ceny,log,nbins,varargin)

  p = inputParser;
  
  p.addRequired('map',@isnumeric);
  p.addRequired('mask',@isnumeric);
  p.addRequired('cenx',@isnumeric);
  p.addRequired('ceny',@isnumeric);
  p.addRequired('log',@isnumeric);
  p.addRequired('nbins',@isnumeric);
  p.addOptional('sig',5,@isnumeric);
  p.addOptional('iter_clip',0,@isnumeric);
  p.addOptional('weight',ones(size(map)),@isnumeric);
  p.addOptional('binedges',[],@isnumeric);
  p.parse(map,mask,cenx,ceny,log,nbins,varargin{:});

  map      = p.Results.map;
  mask     = p.Results.mask;
  cenx     = p.Results.cenx;
  ceny     = p.Results.ceny;
  log      = p.Results.log;
  nbins    = p.Results.nbins;
  sig      = p.Results.sig;
  iter_clip= p.Results.iter_clip;
  weight   = p.Results.weight;
  binedges = p.Results.binedges;
  clear p varargin;
%%

    rad = make_radius_map(map,cenx,ceny);
    if numel(binedges) == 0
        if log
            binedges=logspace(0,log10(max(rad(:))),nbins+1);
            binedges(1)=0;
        else
            binedges=linspace(0,max(rad(:)),nbins+1);
        end
    else
        nbins = numel(binedges) - 1;
    end
  %% setup output parameters

  r=1:nbins;
  prof=r;
  err=r;
  prof2d = rad.*0;
  prof2d_unholy = rad.*0;
  
  npix=r;
  %% loop around the distance map and bin up your signal map

  for i=1:nbins
      spot=(rad >= binedges(i) & rad < binedges(i+1) & mask);
      spotnom=(rad >= binedges(i) & rad < binedges(i+1));
      vecr=rad(spot);
      r(i)=nanmean(vecr(:));
      vecim=map(spot);
      w = weight(spot);
      
      if iter_clip>0
        w = w(vecim==vecim);
        b = vecim(vecim==vecim);
        
        for j=1:iter_clip
            bsp = (b < nanmedian(b(:))+sig*std(b(:)) ... 
                & b > nanmedian(b(:))-sig*std(b(:)));
            b = b(bsp);
            w = w(bsp);
        end
        prof(i)=sum(b.*w)./sum(w);
        
      else
        w = w(vecim==vecim);
        b = vecim(vecim==vecim);
        prof(i)=sum(b.*w)./sum(w);
      end
      
      prof2d(spot) = prof(i);
      prof2d_unholy(spotnom) = prof(i);
      err(i)=nanstd(b)/sqrt(numel(b));
      npix(i)=numel(spot);
   end

   %% get a nice structure to pass through to matlab
  
   profile = struct('r',r,'prof',prof,'err',err,'npix',npix,...
       'prof2d',prof2d,'prof2d_unholy',prof2d_unholy,'binedges',binedges);


return

