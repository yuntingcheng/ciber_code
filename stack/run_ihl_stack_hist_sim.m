function run_ihl_stack_hist_sim(flight,inst,field,hsc_idx,iset,f_ihl,rvir,spire)
% set: 0 - all sim sources in CBmap, 
%      1 - m<22 sources in CBmap

ifield=8;% use SWIRE PSF and masking function

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if strcmp(field,'SIDES')
    if rvir==1
        load(sprintf('%s/stackmapdatsim%d',loaddir,f_ihl*100),'stackmapdatsim');
    else
        load(sprintf('%s/stackmapdatsim%d_rv%d',loaddir,f_ihl*100,rvir),...
            'stackmapdatsim');
    end
elseif strcmp(field,'HSC')
    name = HSC_field_name(hsc_idx);
    if rvir==1
        load(sprintf('%s/stackmapdathsc_%s%d',loaddir,name,f_ihl*100),...
            'stackmapdatsim');
    else
        load(sprintf('%s/stackmapdathsc_%s%d_rv%d',loaddir,name,f_ihl*100,rvir),...
            'stackmapdatsim');
    end
end

dx = 1200;
verbose = 1;
if iset==0
    cbmap = stackmapdatsim(ifield).all.cbmap;
    psmap = stackmapdatsim(ifield).all.psmap;
    mask_inst = stackmapdatsim(ifield).all.mask_inst_clip;
    strmask = stackmapdatsim(ifield).all.strmask;
    strnum = stackmapdatsim(ifield).all.strnum;
elseif iset==1
    cbmap = stackmapdatsim(ifield).sub.cbmap;
    psmap = stackmapdatsim(ifield).sub.psmap;
    mask_inst = stackmapdatsim(ifield).sub.mask_inst_clip;
    strmask = stackmapdatsim(ifield).sub.strmask;
    strnum = stackmapdatsim(ifield).sub.strnum;
end

for im=1:4
    m_min = im + 15;
    m_max = m_min + 1;

%     type = -1;
%     if strcmp(field,'SIDES')
%         [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount]=...
%         stackihl_sim0_hist(flight,inst,ifield,type,m_min,m_max,...
%         dx,cbmap,psmap,mask_inst,strmask,strnum,1,0,verbose,10,[],spire);
%     elseif strcmp(field,'HSC')
%         [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount]=...
%         stackihl_hsc0_hist(flight,inst,ifield,hsc_idx,type,m_min,m_max,...
%         dx,cbmap,psmap,mask_inst,strmask,strnum,1,0,verbose,10,[],spire);
%     end      
%     histdat(im).histcbs = histcb;
%     histdat(im).histpss = histps;
%     histdat(im).Ibinedges_cb = Ibinedges_cb;
%     histdat(im).Ibinedges_ps = Ibinedges_ps;
%     histdat(im).counts = stackcount;
    
    type = 1;
    if strcmp(field,'SIDES')
        [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount]=...
        stackihl_sim0_hist(flight,inst,ifield,type,m_min,m_max,...
        dx,cbmap,psmap,mask_inst,strmask,strnum,1,0,verbose,10,[],spire);
    elseif strcmp(field,'HSC')
        [histcb, histps,Ibinedges_cb,Ibinedges_ps,stackcount]=...
        stackihl_hsc0_hist(flight,inst,ifield,hsc_idx,type,m_min,m_max,...
        dx,cbmap,psmap,mask_inst,strmask,strnum,1,0,verbose,10,[],spire);
    end        
    histdat(im).histcbg = histcb;
    histdat(im).histpsg = histps;
    histdat(im).countg = stackcount;

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst));
if strcmp(field,'SIDES')
    if rvir==1
    save(sprintf('%s/histdatsim_%d',savedir,f_ihl*100),'histdat');
    else
    save(sprintf('%s/histdatsim_%d_rv%d',savedir,f_ihl*100,rvir),...
        'histdat');
    end
elseif strcmp(field,'HSC')
    if rvir==1
    save(sprintf('%s/histdathsc_%s%d',savedir,name,f_ihl*100),'histdat');
    else
    save(sprintf('%s/histdathsc_%s%d_rv%d',savedir,name,f_ihl*100,rvir),...
        'histdat');
    end
end

return