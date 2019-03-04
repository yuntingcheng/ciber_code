function [nmap]=readnoise_realization(ave_n2dps,pixscale,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function take the mean of noise 2DPS of good fields,
%and generate the noise realizatio map using 2 dof chi2 statistics.
%Ref: reduc_adhocnoisemodel.m line 80-101
%
%Input:
%   -ave_n2dps: mean of noise 2DPS
%   -pixscale:beam size (arcsec)
%   -norand(optional):1 turn off the chi2 random
%Output:
%   -nmap: simulated noise map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('ave_n2dps',@isnumeric);
  p.addRequired('pixscale',@isnumeric);
  p.addOptional('norand',0,@isnumeric);
  
  p.parse(ave_n2dps,pixscale,varargin{:});

  ave_n2dps= p.Results.ave_n2dps;
  pixscale = p.Results.pixscale;
  norand = p.Results.norand;
  clear p varargin;

    Cl2din=ave_n2dps.*chi2rnd(ones(1024, 1024).*2)./2;
    %NOTE: divide by 2 b/c 2dof chi2 has mean=2    
    if norand
        Cl2din=ave_n2dps;
    end
    
    nmap=ifft2(fftshift(sqrt(Cl2din)).*fft2(normrnd(0,1,1024)));
    nmap=nmap'./(pixscale/3600.0*pi/180.0);
    nmap=real(nmap)./abs(real(nmap)).*abs(nmap);
end
