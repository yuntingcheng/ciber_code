%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cl_from_Cl2d: This code is nearly same as get_angular_spec.m
%   exept for input is the Cl2d rather than map1 and map2.
%   the only difference is in line 102~109
% Inputs (Required):
%
%   -Cl2d: 2d power spectrum. !!l=0 at the center of map
%   -pixscale: The pixel scale in arcsec.  
%
% Inputs (Optional):
%
%   -nbins:    The number of bins for the resulting power spectrum. 
%                 (Default is 30.)
%   -logscale: 1 for logorithic bins, 0 for linear. (Default is 1.)     
%   -w:        Array of weights for data in fourier space. 
%                 (Default is all ones.) !!!!!in l space- a2d map--YT
%              !!!!note:input weight orientation is not same as Cl2d!!!
%              if we want input weight W(same orientation as Cl2d),              
%              the input should be (fftshift(fftshift(W)))'
% Outputs:
% 
%   -Cl:   The binned power spectrum.
%   -l:    The l-modes correspoding to each bin.
%   -lth:  The angular scale of each bin. (In arcmin.)
%   -dCl:  The standard Knox error bars.
%   -binl: The beginning and ending of each l-bin.
%   -dl:   The size of each l-bin;
%   -Cl2d: The 2d power spectrum.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Cl, l, lth, dCl, binl, dl, Cl2d, lind, nell] = ...
    Cl_from_Cl2d(Cl2d,pixscale,varargin)

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %% Parse data
  p = inputParser;
  
  p.addRequired('Cl2d',@isnumeric);
  p.addRequired('pixscale',@isnumeric);
  p.addOptional('nbins',30,@isnumeric);
  p.addOptional('logbins',1,@isnumeric);
  p.addOptional('w',ones(size(Cl2d)),@isnumeric);
  p.addOptional('superbin',0,@isnumeric);
  p.addOptional('verbose',0,@isnumeric);
  
  p.parse(Cl2d,pixscale,varargin{:});

  Cl2d     = p.Results.Cl2d;
  pixscale = p.Results.pixscale;
  nbins    = p.Results.nbins;
  logbins  = p.Results.logbins;
  w        = p.Results.w;
  superbin = p.Results.superbin;
  verbose  = p.Results.verbose;
  
  clear p varargin;
  %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%
  %%  binning
  
  %% Build multiple frequecy l
  [dim1,dim2] = size(Cl2d);

  l = get_l(dim1,dim2,pixscale);

  %% Make bins either logrithmic or linear.
  if (logbins == 0)
    binl = linspace(0,max(l(:)),nbins);
  else
    binl = logspace(log10(10),log10(max(l(:))),nbins);
  end

  nl = numel(binl);
  
  %%
 %%
  wshft = fftshift(w);
  fftsq=ifftshift(Cl2d').*wshft;
  %%
  
  if verbose == 3
      figure(2); clf; imagesc(Cl2d); caxis([-1e-7,1e-7]); drawnow
      %pause(0.5)
  end
      
  if ~superbin
    nell = zeros(nl-1,1);
    lind = NaN;
      %% Calculate 1d power spectrum
      for i = 1:nl-1
          stuff = fftsq((l >= binl(i)) & (l <= binl(i+1)));
          lower = wshft((l >= binl(i)) & (l <= binl(i+1)));
          if ~isempty(stuff) & isfinite(sum(stuff(lower ~= 0)))    
              Cl(i) = sum(stuff(lower ~= 0))/...
                      sum(lower(lower ~= 0));
              nell(i) = sum(lower(lower ~= 0));
          else
              Cl(i) = 0;
          end
      end
  else 
      
      ellmin = 200;
      whbinl = binl < ellmin;
      whl = find(l < ellmin);
      
      [llo,llosort] = sort(l(whl));
            
      %% Calculate 1d power spectrum
      for i = 1:nl-1
          if binl(i) > ellmin
              stuff = fftsq((l >= binl(i)) & (l <= binl(i+1)));
              lower = wshft((l >= binl(i)) & (l <= binl(i+1)));
              if ~isempty(stuff) & isfinite(sum(stuff(lower ~= 0)))    
                  Cl(i) = sum(stuff(lower ~= 0))/...
                          sum(lower(lower ~= 0));
              else
                  Cl(i) = NaN;
              end
          else
              lowell = binl(i);
              Cl(i) = NaN;
          end
      end
      
      whnonan = isfinite(Cl);
      Cl = Cl(whnonan);
      binl = [lowell,binl(whnonan)];

      Cllo = fftsq(whl(llosort)) ./ wshft(whl(llosort));
      lind = whl(llosort);
      
      %Cl = [Cllo(2:end)',Cl];
      %binl = [llo',binl];
  end

  %% Return cosmic variance error.
  %% The formula for the cosmic variance error can be expressed in terms of
  %% the multipole moment \delta C_l = \sqrt(2/(2*l+1)/(\delta l * f_sky)) where
  %% fsky = (square degrees in map)/(square degrees in sky) and the total
  %% square degrees in the sky is 4 \pi * (180/4\pi)^2.
  sqdeg_in_sky = 4*pi*(180/pi)^2;
  sqdeg_in_map = (dim1*pixscale/3600*dim2*pixscale/3600);
  fsky = sqdeg_in_map/sqdeg_in_sky;

  l = (binl(2:end)+binl(1:end-1))/2.0;
  bin = (binl(2:end)-binl(1:end-1));

  if superbin
      
      Cl = [Cllo(1),(Cllo(2) + Cllo(5))./2,(Cllo(3)+Cllo(4))./2,Cl];
      l = [llo(1:3)',l];
      bin = diff([l,2.*l(end)]);
      lind = lind(1:3);

  end
  
  dCl = Cl.*sqrt(2.0./(2.0*l+1.0)./bin/fsky);
  
  %% Arcmin.
  lth = 2*60*180.0./l;
  dl = bin;

return

