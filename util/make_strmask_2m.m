function [strmask,num]=make_strmask_2m(flight,inst,ifield,...
                          alpha,beta,m_max,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce a star mask from catalog
% the masking radius is defined by r= alpha *m + beta
%
% the magnitude data is in AB magnitude, and linear interpolate to 
% CIBER I and H effective wavelength. 
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - alpha: radius slope [arcsec/mag]
% - beta: radius intersection [arcsec]
% - m_max: max masking magnitude
%(Optional)
% - catname: catalog name - 'PSC','XSC','PSCrej','XSCrej', default:'PSC'
% - rel: reliability in 2M rej cat - 'A'-'F', default:'A'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('alpha',@isnumeric);
  p.addRequired('beta',@isnumeric);
  p.addRequired('m_max',@isnumeric);
  p.addOptional('rel','A',@isstr);
  p.addOptional('catname','PSC',@isstr);
  
  p.parse(flight,inst,ifield,alpha,beta,m_max,varargin{:});

  flight     = p.Results.flight;
  inst     = p.Results.inst;
  ifield = p.Results.ifield;
  alpha    = p.Results.alpha;
  beta  = p.Results.beta;
  m_max        = p.Results.m_max;
  catname = p.Results.catname;
  rel  = p.Results.rel;
  
  clear p varargin;

catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/',catname,'/');

quad_arr=['A','B','C','D'];
dt=get_dark_times(flight,inst,ifield);

strmask=zeros(1024);
numk=zeros(1024);
for iquad=1:4

%%% read cat data %%%
quad=quad_arr(iquad);
if strcmp(catname,'PSCrej')
    catfile=strcat(catdir,dt.name,'_',quad,'_2m_',rel,'.txt');
else
    catfile=strcat(catdir,dt.name,'_',quad,'_2m.txt');
end

M = csvread(catfile,1);

x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');

x_arr=x_arr+1;
y_arr=y_arr+1;

if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end

if strcmp(catname,'XSC') || strcmp(catname,'XSCrej') 
    rext_arr=squeeze(M(:,11)'); %[arcsec]
    rext_arr=rext_arr./7;% [pix]
end
  
%%% produce mask %%%%
%rad_arr=alpha.*m_arr+beta; % [arcsec]

%rad_arr=1930.*m_arr.^-0.8053-172.6; % [arcsec]
rad_arr = -41.6.*exp(-(m_arr-7.76).^2./(1.298).^2) + ...
    1455.*exp(-(m_arr+12.82).^2./(13.94).^2); % [arcsec]

rad_arr=rad_arr./7; % [pix]


% in case r goes to negative before m_max
rad_arr(find(rad_arr<1))=1;

% if mask radius < rext from 2M cat, set radius to r_ext
if strcmp(catname,'XSC') || strcmp(catname,'XSCrej') 
    sp=find(rad_arr<rext_arr);
    rad_arr(sp)=rext_arr(sp);
end

sp=find(m_arr<m_max);
subm_arr=m_arr(sp);
subrad_arr=rad_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

mask=ones(512);
qnum=zeros(512);

for i=1:numel(subm_arr)
    radmap=make_radius_map(mask,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < subrad_arr(i));
    mask(sp1)=0;
    qnum(sp1)=qnum(sp1)+1;
end


%%% stick to large mask %%%%
if iquad==1
    strmask(1:512,1:512)=mask;
    num(1:512,1:512)=qnum;
elseif iquad==2
    strmask(513:1024,1:512)=mask;
    num(513:1024,1:512)=qnum;
elseif iquad==3
    strmask(1:512,513:1024)=mask;
    num(1:512,513:1024)=qnum;
elseif iquad==4
    strmask(513:1024,513:1024)=mask;
    num(513:1024,513:1024)=qnum;
end

end

end
