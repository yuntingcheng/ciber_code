function stack_preprocess_hsc(flight,inst)

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdathsc'),'maskdat');
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

for hsc_idx=0:11
    %%
    [name,mask_inst] = HSC_fields_info(hsc_idx);

    map1 = fits_read(strcat(srcmapdir,'srcmap_sim1_all_hsc_',name,'.fits'));
    map2 = fits_read(strcat(srcmapdir,'srcmap_sim2_all_hsc_',name,'.fits'));
    
    %{
    figure
    setwinsize(gcf,1200,300)
    subplot(1,3,1)
    imageclip(map1);
    title('m<20');
    h = colorbar;
    ylabel(h, 'nW/m^2/sr');
    subplot(1,3,2)
    imageclip(map2);
    h = colorbar;
    ylabel(h, 'nW/m^2/sr');
    title('m>20');
    subplot(1,3,3)
    imageclip(map2.*mask_inst);
    title('m>20 masked');
    h = colorbar;
    ylabel(h, 'nW/m^2/sr');
    %}
    
    psmap_raw = map1;
    cbmap_raw = map1 + map2;
    
    strmask = maskdat.mask(hsc_idx+1).m_max(20).strmask_stack;
    strnum = maskdat.mask(hsc_idx+1).m_max(20).strnum_stack;
    totmask = mask_inst.*strmask;
    
    sigmask = totmask;
    cbmean = mean(cbmap_raw(find(sigmask)));
    psmean = mean(psmap_raw(find(sigmask)));
    
    p = polyfitweighted2(1:1024,1:1024,cbmap_raw,2,sigmask);
    polymapcb = polyval2(p,1:1024,1:1024);
    cbmap = cbmap_raw - polymapcb;
    cb_bk = mean(cbmap(find(sigmask)));
    cbmap = cbmap - cb_bk;

    p = polyfitweighted2(1:1024,1:1024,psmap_raw,2,sigmask);
    polymapps = polyval2(p,1:1024,1:1024);
    psmap = psmap_raw - polymapps;
    ps_bk = mean(psmap(find(sigmask)));
    psmap = psmap - ps_bk;

    %%% write the data %%%
    stackmapdat(hsc_idx+1).cbmap = cbmap;
    stackmapdat(hsc_idx+1).psmap = psmap;
    stackmapdat(hsc_idx+1).mask_inst_clip = mask_inst;
    stackmapdat(hsc_idx+1).strmask = strmask;
    stackmapdat(hsc_idx+1).strnum = strnum;
    stackmapdat(hsc_idx+1).m_min_arr = 16:19;
    stackmapdat(hsc_idx+1).m_max_arr = 17:20;
    stackmapdat(hsc_idx+1).cbmean = cbmean;
    stackmapdat(hsc_idx+1).psmean = psmean;
    figure
    setwinsize(gcf,1000,300)
    subplot(1,2,1)
    histogram(cbmap.*mask_inst.*strmask);
    subplot(1,2,2)
    histogram(psmap.*mask_inst.*strmask);
    title(hsc_idx);
    %%% get mask_inst_clip from strmask with limits %%%
    for m_max=17:20
        strmask = maskdat.mask(hsc_idx+1).m_max(m_max).strmask_stack;
        strnum = maskdat.mask(hsc_idx+1).m_max(m_max).strnum_stack;
        stackmapdat(hsc_idx+1).m_max(m_max).strmask = strmask;
        stackmapdat(hsc_idx+1).m_max(m_max).strnum = strnum;
    end
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackmapdathsc',savedir),...
    'stackmapdat');
return