function run_ihl_stack_hist(flight,inst,ifield,spire,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('spire',@islogical);
  p.addOptional('rmin',nan,@isnumeric);
  
  p.parse(flight,inst,ifield,spire,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  spire    = p.Results.spire;
  rmin     = p.Results.rmin;
  
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if rmin==2
    load(sprintf('%sstackmapdat_rmin2',loaddir),'stackmapdat');
elseif isnan(rmin)
    load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
end
dx = 1200;
verbose = 1;
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
return_count=false;

if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end

for im= 10:13
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    type = -1;
    [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount,~,~,~]=...
    stackihl_ps0_hist(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,1,0,verbose,10,stackband,[],spire,rmin,return_count);
    histdat(im).histcbs = histcb;
    histdat(im).histpss = histps;
    histdat(im).Ibinedges_cb = Ibinedges_cb;
    histdat(im).Ibinedges_ps = Ibinedges_ps;
    histdat(im).counts = stackcount;
    
    type = 1;
    [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount,~,~,~]=...
    stackihl_ps0_hist(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,1,0,verbose,10,stackband,[],spire,rmin,return_count);
    histdat(im).histcbg = histcb;
    histdat(im).histpsg = histps;
    histdat(im).countg = stackcount;

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

if rmin==2
    save(sprintf('%s/histdat_%s_rmin2',savedir,dt.name),'histdat');
elseif isnan(rmin)
    save(sprintf('%s/histdat_%s',savedir,dt.name),'histdat');
end

return