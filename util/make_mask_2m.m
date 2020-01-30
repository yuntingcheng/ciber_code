function [mask,num]=make_mask_2m(flight,inst,ifield,m_min,m_max,Ith,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the mask from 2MASS
%Use mag extrapolated to PanSTARRS y band 
%in order to be consistent with PS mask
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8
% - band: j,h,k(2M),I,H(CB-interpolated),y(PS-interpolated)
% - m_min: min masking magnitude (PS y band extrapolated)
% - m_max: max masking magnitude (PS y band extrapolated)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('m_min',@isnumeric);
  p.addRequired('m_max',@isnumeric);
  p.addRequired('Ith',@isnumeric);
  p.addOptional('rmin',nan,@isnumeric);
  p.addOptional('PSmatch',nan,@isnumeric);
  p.addOptional('verbose',true,@islogical);
  p.parse(flight,inst,ifield,m_min,m_max,Ith,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  m_min    = p.Results.m_min;
  m_max    = p.Results.m_max;
  Ith      = p.Results.Ith;
  rmin     = p.Results.rmin;
  PSmatch  = p.Results.PSmatch;
  verbose  = p.Results.verbose;
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/PSC/');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.csv');

M = csvread(catfile,1);

if inst==1
    x_arr=squeeze(M(:,4)');
    y_arr=squeeze(M(:,3)');
else
    x_arr=squeeze(M(:,6)');
    y_arr=squeeze(M(:,5)');
end

x_arr=x_arr+1;
y_arr=y_arr+1;

% keySet =   {'j','h','k','I','H','y'};
% idxSet = [5,6,7,8,9,10];
% idxObj = containers.Map(keySet,idxSet);
% m_arr=squeeze(M(:,idxObj(band))');
if inst==1
    m_arr=squeeze(M(:,10)');
else
    m_arr=squeeze(M(:,11)');
end

mI_arr = squeeze(M(:,10)');
if PSmatch~=0
    match_arr = squeeze(M(:,13)');
    sp=find(mI_arr<=m_max & mI_arr>m_min & match_arr==0);
else
    sp=find(mI_arr<=m_max & mI_arr>m_min);
end

subm_arr=m_arr(sp);
submI_arr=mI_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

if isnan(rmin)
%     rad_arr = get_mask_radius_th(inst,ifield,subm_arr,Ith);
    rad_arr = get_mask_radius_th(1,ifield,submI_arr,Ith);
else
%     rad_arr = get_mask_radius_th(inst,ifield,subm_arr,Ith,'rmin',rmin);
    rad_arr = get_mask_radius_th(1,ifield,submI_arr,Ith,'rmin',rmin);
end

mask = ones(1024);
num = zeros(1024);

if numel(subm_arr)>100
    idx_print=floor(numel(subm_arr)/100) * (1:100);
else
    idx_print=[];
end
print_count=0;

for i=1:numel(subm_arr)
    radmap=make_radius_map(mask,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7);
    mask(sp1)=0;
    num(sp1) = num(sp1) + 1;
    
    if verbose
        if ismember(i,idx_print)
            print_count = print_count + 1;
            disp(sprintf('making mask for %d %% sources',print_count));
        end
    end

end

end
