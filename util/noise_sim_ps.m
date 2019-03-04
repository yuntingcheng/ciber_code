function nsimps=noise_sim_ps(flight,inst,field,mask,mkk,weight,diff,g1,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do photon noise + shot noise simulation many times and save the 1DPS of
% all the simulation 
%Input:
%   -flgiht:40030/36277/36265
%   -inst:1/2
%   -field:field number (see get_fields.m)
%   -mask:mask(use flight mask)
%   -mkk: mkk of that mask
%   -diff: diff(1) or full integration(0)
%   -g1: cal factor g1, NOTE:g1 must be negative
%   -pixscale(optional):pixel scale in arcsec (default:7)
%   -weight(optional):Fourier weight(default:ones(1024))
%   -nsim(optional):# of simulation, also size of output PS (default:100)
%Output: 
%   -nsimps.l:l
%   -nsimps.rnCl_arr: read noise PS (size:nsim, numel(l))
%   -nsimps.phCl_arr: photon noise PS
%   -nsimps.nCl_arr: noise PS = rnCl_arr+phCl_arr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Parse data
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('field',@isnumeric);
  p.addRequired('mask',@isnumeric);
  p.addRequired('mkk',@isnumeric);
  p.addRequired('weight',@isnumeric);
  p.addRequired('diff',@isnumeric);  
  p.addRequired('g1',@isnumeric);  
  p.addOptional('pixscale',7,@isnumeric);
  p.addOptional('nsim',100,@isnumeric);
  
  p.parse(flight,inst,field,mask,mkk,weight,diff,g1,varargin{:});

  flight     = p.Results.flight;
  inst     = p.Results.inst;
  field     = p.Results.field;
  mask     = p.Results.mask;
  mkk     = p.Results.mkk;
  weight     = p.Results.weight;
  diff     = p.Results.diff;
  g1     = p.Results.g1;
  pixscale     = p.Results.pixscale;
  nsim     = p.Results.nsim;
  
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,~,l]=get_Cl(randn(size(mask)),mask,mkk,pixscale,weight);

%%% get Fourier weight by stack Cl2d
nCl2d_arr=zeros(nsim,1024,1024);
nmap_arr=zeros(nsim,1024,1024);
rnmap_arr=zeros(nsim,1024,1024);
shotmap_arr=zeros(nsim,1024,1024);
for j=1:nsim
[nmap,rnmap,shotmap]=noise_realization(flight,inst,field,diff,g1);
[~,~,nCl2d]=get_Cl(nmap,mask,mkk,pixscale,weight);
nCl2d_arr(j,:,:)=nCl2d;
nmap_arr(j,:,:)=nmap;
rnmap_arr(j,:,:)=rnmap;
shotmap_arr(j,:,:)=shotmap;

fprintf('field=%d,jsim=%d\n',field,j);
end
meanCl2d=squeeze(mean(nCl2d_arr));
stdCl2d=squeeze(std(nCl2d_arr));
weight_sim=(fftshift(fftshift(1./stdCl2d)))';
clear nCl2d_arr


%%% now get 1DPS with weight based on sim 2DPS
phClws_arr=zeros(nsim,numel(l));
rnClws_arr=zeros(nsim,numel(l));
nClws_arr=zeros(nsim,numel(l));
for j=1:nsim
nmap=squeeze(nmap_arr(j,:,:));
rnmap=squeeze(rnmap_arr(j,:,:));
shotmap=squeeze(shotmap_arr(j,:,:));
nCl=get_Cl(nmap,mask,mkk,pixscale,weight_sim);
rnCl=get_Cl(rnmap,mask,mkk,pixscale,weight_sim);
phCl=get_Cl(shotmap,mask,mkk,pixscale,weight_sim);
nClws_arr(j,:)=nCl;rnClws_arr(j,:)=rnCl;phClws_arr(j,:)=phCl;
end

%%% now get 1DPS with input weight
phCl_arr=zeros(nsim,numel(l));
rnCl_arr=zeros(nsim,numel(l));
nCl_arr=zeros(nsim,numel(l));
for j=1:nsim
nmap=squeeze(nmap_arr(j,:,:));
rnmap=squeeze(rnmap_arr(j,:,:));
shotmap=squeeze(shotmap_arr(j,:,:));
nCl=get_Cl(nmap,mask,mkk,pixscale,weight);
rnCl=get_Cl(rnmap,mask,mkk,pixscale,weight);
phCl=get_Cl(shotmap,mask,mkk,pixscale,weight);
nCl_arr(j,:)=nCl;rnCl_arr(j,:)=rnCl;phCl_arr(j,:)=phCl;
end
clear nmap_arr rnmap_arr shotmap_arr



nsimps.l=l;
nsimps.rnCl_arr=rnCl_arr;
nsimps.phCl_arr=phCl_arr;
nsimps.nCl_arr=nCl_arr;
nsimps.rnClws_arr=rnClws_arr;
nsimps.phClws_arr=phClws_arr;
nsimps.nClws_arr=nClws_arr;
nsimps.meanCl2d=meanCl2d;
nsimps.stdCl2d=stdCl2d;
nsimps.weight_sim=weight_sim;
return