function run_ihl_stack_map_hsc(flight,inst,hsc_idx,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('hsc_idx',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  p.addOptional('sample_type','jack_random',@ischar);
  p.addOptional('subpix',true,@islogical);
  
  p.parse(flight,inst,hsc_idx,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  hsc_idx   = p.Results.hsc_idx;
  masklim = p.Results.masklim;
  sample_type=p.Results.sample_type;
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
clipmaxs = ones([2,nbins])*inf;
clipmins = -ones([2,nbins])*inf;
if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end
%%
for im= 1:numel(m_min_arr)

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    if masklim
        mask_inst = stackmapdat(hsc_idx+1).mask_inst_clip;
        strmask = stackmapdat(hsc_idx+1).m_max(m_max).strmask;
        strnum = stackmapdat(hsc_idx+1).m_max(m_max).strnum;
    else
        mask_inst = stackmapdat(hsc_idx+1).mask_inst_clip;
        strmask = stackmapdat(hsc_idx+1).strmask;
        strnum = stackmapdat(hsc_idx+1).strnum;    
    end
    
    stackdat.m_min = m_min;
    stackdat.m_max = m_max;
    
    srcdat = ps_src_select(flight,inst,hsc_idx,m_min,m_max,mask_inst,stackband,...
    'sample_type',sample_type,'Nsub',10,'HSC',true);

    stackdat.r_arr = rbins;
    for isub=1:10
        [~,~,~,profcbg,profpsg,profhitg] = ...
            stackihl_ps0_hist_map(flight,inst,8,dx,cbmap,psmap,...
            mask_inst,strmask,strnum,1,verbose,nan,clipmaxs,...
            clipmins,srcdat.sub(isub).xg_arr,...
            srcdat.sub(isub).yg_arr,srcdat.sub(isub).mg_arr,subpix);

        fprintf('stack %s, %d<m<%d, isub %d\n',...
            name,m_min,m_max,isub);
        
        profcbg(profhitg==0) = 0;
        profpsg(profhitg==0) = 0;
        stackdat.sub(isub).countg = srcdat.sub(isub).Ng;
        stackdat.sub(isub).profcbg = profcbg;
        stackdat.sub(isub).profpsg = profpsg;
        stackdat.sub(isub).profhitg = profhitg;
    end

    sp100 = find(stackdat.r_arr>100);

    %%% profile combining all subset
    profcbg = zeros(size(stackdat.r_arr));
    profpsg = zeros(size(stackdat.r_arr));
    profhitg = zeros(size(stackdat.r_arr));
    countg = 0;
    for isub=1:10
        profcbg = profcbg + ...
            stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        profpsg = profpsg + ...
            stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        profhitg = profhitg + stackdat.sub(isub).profhitg;
        countg = countg + stackdat.sub(isub).countg;
    end
    stackdat.all.profcbg = profcbg./profhitg;
    stackdat.all.profpsg = profpsg./profhitg;
    stackdat.all.countg = countg;
    stackdat.all.profcbg100 = sum(profcbg(sp100))./sum(profhitg(sp100));
    stackdat.all.profpsg100 = sum(profpsg(sp100))./sum(profhitg(sp100));

    %%% profile of jackknife samples (leave one out)
    for isub=1:10
        jackcbg = profcbg - ...
            stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        jackpsg = profpsg - ...
            stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        jackhitg = profhitg - stackdat.sub(isub).profhitg;
        stackdat.jack(isub).profcbg = jackcbg./jackhitg;
        stackdat.jack(isub).profpsg = jackpsg./jackhitg;
        stackdat.jack(isub).profcbg100 = ...
            sum(jackcbg(sp100))./sum(jackhitg); 
        stackdat.jack(isub).profpsg100 = ...
            sum(jackpsg(sp100))./sum(jackhitg); 
    end

    %%% error bar with jackknife
    errcbg = zeros(size(stackdat.r_arr));
    errpsg = zeros(size(stackdat.r_arr));
    errcbg100 = 0;
    errpsg100 = 0;
    for isub=1:10
    errcbg = errcbg + ...
        (stackdat.jack(isub).profcbg - stackdat.all.profcbg).^2;
    errpsg = errpsg + ...
        (stackdat.jack(isub).profpsg - stackdat.all.profpsg).^2;
    errcbg100 = errcbg100 + ...
        (stackdat.jack(isub).profcbg100 - stackdat.all.profcbg100).^2;
    errpsg100 = errpsg100 + ...
        (stackdat.jack(isub).profpsg100 - stackdat.all.profpsg100).^2;
    end
    stackdat.errjack.profcbg = sqrt(errcbg.*(9/10));
    stackdat.errjack.profpsg = sqrt(errpsg.*(9/10));
    stackdat.errjack.profcbg100 = sqrt(errcbg100.*(9/10));
    stackdat.errjack.profpsg100 = sqrt(errpsg100.*(9/10));
    
    stackdathsc(im).stackdat = stackdat;
    clear stackdat

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if masklim
    save(sprintf('%s/stackdathsc_%s_masklim',...
        savedir,name),'stackdathsc');        
else
    save(sprintf('%s/stackdathsc_%s',...
        savedir,name),'stackdathsc');
end
clear stackdat

    
return
