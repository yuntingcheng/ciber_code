function run_ihl_stack_bk_hist(flight,inst,ifield_arr,run,spire,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield_arr',@isnumeric);
  p.addRequired('run',@isnumeric);
  p.addRequired('spire',@islogical);
  p.addOptional('rmin',nan,@isnumeric);
  p.addOptional('load_count',false,@isnumeric);
  
  p.parse(flight,inst,ifield_arr,run,spire,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  run   = p.Results.run;
  spire    = p.Results.spire;
  rmin     = p.Results.rmin;
  load_count= p.Results.load_count;
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if rmin==2
    load(sprintf('%sstackmapdat_rmin2',loaddir),'stackmapdat');
elseif isnan(rmin)
    load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
end
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

dx = 1200;
for ifield=ifield_arr
    dt=get_dark_times(flight,inst,ifield);    
    cbmap = stackmapdat(ifield).cbmap;
    psmap = stackmapdat(ifield).psmap;
    mask_inst = stackmapdat(ifield).mask_inst_clip;
    strmask = stackmapdat(ifield).strmask;
    
    counts_arr = [];
    countg_arr = [];
    if load_count
        if rmin==2
            load(sprintf('%shistdat_%s_rmin2',loaddir,dt.name),'histdat');
        elseif isnan(rmin)
            load(sprintf('%shistdat_%s',loaddir,dt.name),'histdat');
        end
        for im = 10:13
            counts_arr = [counts_arr, histdat(im).counts];
            countg_arr = [countg_arr, histdat(im).countg];
        end
        Ibinedges_cb = histdat(im).Ibinedges_cb;
        Ibinedges_ps = histdat(im).Ibinedges_ps;
    else
        m_min_arr = stackmapdat(ifield).m_min_arr;
        m_max_arr = stackmapdat(ifield).m_max_arr;
        strnum = stackmapdat(ifield).strnum;
        if inst == 1
            stackband = 'I';
        else
            stackband = 'H';
        end
        for im=10:13
            m_min = m_min_arr(im);
            m_max = m_max_arr(im);
            [~,~,Ibinedges_cb,Ibinedges_ps,stackcount,~,~,~]=...
            stackihl_ps0_hist(flight,inst,ifield,-1,m_min,m_max,dx,cbmap,psmap,...
            mask_inst,strmask,strnum,1,0,false,10,stackband,[],spire,rmin,true);
            counts_arr = [counts_arr, stackcount];
            [~,~,Ibinedges_cb,Ibinedges_ps,stackcount,~,~,~]=...
            stackihl_ps0_hist(flight,inst,ifield,1,m_min,m_max,dx,cbmap,psmap,...
            mask_inst,strmask,strnum,1,0,false,10,stackband,[],spire,rmin,true);
            countg_arr = [countg_arr, stackcount];
        end
    end
    N_arr = unique([counts_arr, countg_arr]);
    N_arr = N_arr(find(N_arr>0));
    for iter=(0:4)*10 + run
        [histbkdat]=stackihl_ps_randomN_hist(N_arr,dx,cbmap,psmap,...
            Ibinedges_cb,Ibinedges_ps,mask_inst.*strmask,1024,spire);
        
        if rmin==2
            save(strcat(savedir,'bk_ps/',dt.name,...
                '_histbkdat','_',num2str(iter),'_rmin2'),'histbkdat');
        elseif isnan(rmin)
            save(strcat(savedir,'bk_ps/',dt.name,...
                '_histbkdat','_',num2str(iter)),'histbkdat');
        end      
    end
end

return