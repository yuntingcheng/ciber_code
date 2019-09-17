function run_ihl_stack_map_hsc_bk(flight,inst,hsc_idx,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('hsc_idx',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  p.addOptional('subpix',false,@islogical);
  
  p.parse(flight,inst,hsc_idx,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  hsc_idx  = p.Results.hsc_idx;
  masklim = p.Results.masklim;
  subpix   = p.Results.subpix;
  clear p varargin;
%%
mypaths=get_paths(flight);
[name,~] = HSC_fields_info(hsc_idx);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdathsc',loaddir),'stackmapdat');
dx = 1200;
verbose = false;
cbmap = stackmapdat(hsc_idx+1).cbmap;
psmap = stackmapdat(hsc_idx+1).psmap;
m_min_arr = stackmapdat(hsc_idx+1).m_min_arr;
m_max_arr = stackmapdat(hsc_idx+1).m_max_arr;

nbins = 25;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges = profile.binedges;
rbins = binedges2bins(rbinedges).*0.7;

N_arr = zeros([1,numel(m_min_arr)]);
for im= 1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    mask_inst = stackmapdat(hsc_idx+1).mask_inst_clip;
    srcdat = ps_src_select(flight,inst,hsc_idx,m_min,m_max,...
        mask_inst,'HSC',true);
    N_arr(im) = srcdat.Ng;
end

for iter=1:100
    stackdat.r_arr = rbins;
    for im= 1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
        stackdat(im).m_min=m_min;
        stackdat(im).m_max=m_max;

        mask_inst = stackmapdat(hsc_idx+1).mask_inst_clip;
        if masklim
            strmask = stackmapdat(hsc_idx+1).m_max(m_max).strmask;
        else
            strmask = stackmapdat(hsc_idx+1).strmask;
        end

        [profcb_arr,profps_arr,hitmap_arr]=stackihl_ps0_hist_map_bk...
        (dx,cbmap,psmap,mask_inst,strmask,N_arr(im),verbose,subpix);

        profcb_arr(hitmap_arr==0) = 0;
        profps_arr(hitmap_arr==0) = 0;

        stackdat(im).Ng = N_arr(im);
        stackdat(im).profcbg_arr = profcb_arr;
        stackdat(im).profpsg_arr = profps_arr;
        stackdat(im).hitmapg_arr = hitmap_arr;
        fprintf('stack %s, %d<m<%d, %d src, iter %d\n',...
            name,m_min,m_max,N_arr(im),iter);
    end

        stackdathscbk(iter).stackdat = stackdat;
        clear stackdat
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if masklim
    save(sprintf('%s/hsc/stackdathsc_%s_masklim_bk',...
        savedir,name),'stackdathscbk');  
else
    save(sprintf('%s/hsc/stackdathsc_%s_bk',...
        savedir,name),'stackdathscbk');
end
clear stackhscbk

return
