function run_ihl_stack_map_bk(flight,inst,ifield,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addOptional('subpix',false,@isnumeric);
  p.addOptional('masklim',false,@islogical);
  
  p.parse(flight,inst,ifield,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  subpix   = p.Results.subpix;
  masklim = p.Results.masklim;
  clear p varargin;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end

load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dx = 1200;
verbose = false;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
m_min_arr = m_min_arr(10:13);
m_max_arr = m_max_arr(10:13);


savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/bk_ps/'));
%%
nbins = 25;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges = profile.binedges;
rbins = binedges2bins(rbinedges).*0.7;

Ng_arr = zeros([1,numel(m_min_arr)]);
Ns_arr = zeros([1,numel(m_min_arr)]);
for im= 1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    if masklim
        mask_inst = stackmapdat(ifield).m_max(m_max).mask_inst_clip;
    else
        mask_inst = stackmapdat(ifield).mask_inst_clip;
    end
    srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,stackband);
    Ng_arr(im) = srcdat.Ng;
    Ns_arr(im) = srcdat.Ns;
end


for iter=1:100
stackdat.r_arr = rbins;
    for im= 1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
        stackdat(im).m_min=m_min;
        stackdat(im).m_max=m_max;

        if masklim
            mask_inst = stackmapdat(ifield).m_max(m_max).mask_inst_clip;
            strmask = stackmapdat(ifield).m_max(m_max).strmask;
        else
            mask_inst = stackmapdat(ifield).mask_inst_clip;
            strmask = stackmapdat(ifield).strmask;
        end

        [profcb_arr,profps_arr,hitmap_arr]=stackihl_ps0_hist_map_bk...
        (dx,cbmap,psmap,mask_inst,strmask,[Ns_arr(im),Ng_arr(im)],...
        verbose,subpix);
        profcb_arr(hitmap_arr==0) = 0;
        profps_arr(hitmap_arr==0) = 0;

        stackdat(im).Ns = Ns_arr(im);
        stackdat(im).Ng = Ng_arr(im);
        stackdat(im).profcbs_arr = profcb_arr(1,:);
        stackdat(im).profpss_arr = profps_arr(1,:);
        stackdat(im).hitmaps_arr = hitmap_arr(1,:);
        stackdat(im).profcbg_arr = profcb_arr(2,:);
        stackdat(im).profpsg_arr = profps_arr(2,:);
        stackdat(im).hitmapg_arr = hitmap_arr(2,:);
        fprintf('stack %s, %d<m<%d, %d stars %d gals, iter %d\n',...
        dt.name,m_min,m_max,Ns_arr(im),Ng_arr(im),iter);
    end

    stackdatbk(iter).stackdat = stackdat;
    clear stackdat
end

if masklim
    save(sprintf('%s/stackdat_%s_masklim_bk',...
        savedir,dt.name),'stackdatbk');  
else
    save(sprintf('%s/stackdat_%s_bk',...
        savedir,dt.name),'stackdatbk');
end
clear stackdatbk

%%
return
