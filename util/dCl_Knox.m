function [dCl]=dCl_Knox(Cl,binl,beam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given Cl and binl from get_angular_spec, calculate Knox err dCl.
%Adopt from get_angular_spec.m
%Input:
%   -Cl
%   -binl
%   -beam
%Output:
%   -dCl: same as dCl from get_angular_spec.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sqdeg_in_sky = 4*pi*(180/pi)^2;
  sqdeg_in_map = (1024*beam/3600*1024*beam/3600);
  fsky = sqdeg_in_map/sqdeg_in_sky;

  l = (binl(2:end)+binl(1:end-1))/2.0;
  bin = (binl(2:end)-binl(1:end-1));

  l=l(1:numel(Cl));bin=bin(1:numel(Cl));
  dCl = Cl.*sqrt(2.0./(2.0*l+1.0)./bin/fsky);
end
