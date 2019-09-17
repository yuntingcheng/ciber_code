function [srcdat]=tm_src_select(flight,inst,ifield,m_min,m_max,...
    mask_inst,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select sources to be stack
% sample_type:
% - 'all': all the sources
% - 'jack_random': randomly separate sources into Nsub subsets
% - 'jack_region': separate sources into 16(4x4) spatial regions
% - 'jack_random_uni': randomly separate uniformizec sources into Nsub subsets
% - 'jack_region_uni': separate uniformized sources into 16(4x4) spatial regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('m_min',@isnumeric);
  p.addRequired('m_max',@isnumeric);
  p.addRequired('mask_inst',@isnumeric);
  p.addOptional('Nsub',16,@isnumeric);
  p.addOptional('sample_type','all',@ischar);
  p.parse(flight,inst,ifield,m_min,m_max,mask_inst,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  m_min    = p.Results.m_min;
  m_max    = p.Results.m_max;
  mask_inst= p.Results.mask_inst;
  sample_type = p.Results.sample_type;
  Nsub     = p.Results.Nsub;
  clear p varargin;
%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/PSC/');
dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.csv');

M = csvread(catfile,1);
x1_arr=squeeze(M(:,4)');
y1_arr=squeeze(M(:,3)');
x2_arr=squeeze(M(:,6)');
y2_arr=squeeze(M(:,5)');    

x1_arr=x1_arr+1;
y1_arr=y1_arr+1;
x2_arr=x2_arr+1;
y2_arr=y2_arr+1;

m_arr=squeeze(M(:,7)');% select by 2MASS j band
    
sp=find(x1_arr>0.5 & x1_arr<1024.5 & y1_arr>0.5 & y1_arr<1024.5 & ...
    x2_arr>0.5 & x2_arr<1024.5 & y2_arr>0.5 & y2_arr<1024.5);

x1_arr = x1_arr(sp);
y1_arr = y1_arr(sp);
x2_arr = x2_arr(sp);
y2_arr = y2_arr(sp);
m_arr = m_arr(sp);

%%% count the center pix map
sp=find(m_arr<=m_max & m_arr>m_min);

m_arr=m_arr(sp);
x1_arr=x1_arr(sp);
y1_arr=y1_arr(sp);
x2_arr=x2_arr(sp);
y2_arr=y2_arr(sp);

%%% select sources not coexist with others in the same pixel %%%
mask_inst1 = squeeze(mask_inst(1,:,:));
mask_inst2 = squeeze(mask_inst(2,:,:));
subm_arr=[];
subx1_arr=[];
suby1_arr=[];
subx2_arr=[];
suby2_arr=[];
for i=1:numel(sp)
    if  mask_inst1(round(x1_arr(i)),round(y1_arr(i)))==1 ...
            & mask_inst2(round(x2_arr(i)),round(y2_arr(i)))==1
        subm_arr=[subm_arr,m_arr(i)];
        subx1_arr=[subx1_arr,x1_arr(i)];
        suby1_arr=[suby1_arr,y1_arr(i)];
        subx2_arr=[subx2_arr,x2_arr(i)];
        suby2_arr=[suby2_arr,y2_arr(i)];
    end
end

randidx = randperm(numel(subm_arr));
if inst==1
    xs_arr = subx1_arr(randidx);
    ys_arr = suby1_arr(randidx);
else
    xs_arr = subx2_arr(randidx);
    ys_arr = suby2_arr(randidx);        
end
    ms_arr = subm_arr(randidx);
%%
if strcmp(sample_type,'all')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ns = numel(xs_arr);
    srcdat.xs_arr = xs_arr;
    srcdat.ys_arr = ys_arr;
    srcdat.ms_arr = ms_arr;
    return
end
%%
if strcmp(sample_type,'jack_random')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ns = numel(xs_arr);
    for isub=1:Nsub
        srcdat.sub(isub).xs_arr = xs_arr(isub:Nsub:numel(xs_arr));
        srcdat.sub(isub).ys_arr = ys_arr(isub:Nsub:numel(xs_arr));
        srcdat.sub(isub).ms_arr = ms_arr(isub:Nsub:numel(xs_arr));
        srcdat.sub(isub).Ns = numel(srcdat.sub(isub).xs_arr);
    end
    return
end
%%

return
