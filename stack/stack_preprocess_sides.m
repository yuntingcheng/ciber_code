function stack_preprocess_sides(flight,inst,ifield)

mypaths=get_paths(flight);

dt=get_dark_times(flight,inst,ifield);

%%%%%%%%%%%% srcmaps %%%%%%%%%%%%%%%
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

cbmap_raw = zeros(720);
for m_min = 13:28
    m_max = m_min + 1;
    map = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
        num2str(m_min),'_',num2str(m_max),'sides.fits'));
    cbmap_raw = map + cbmap_raw;
end

psmap_raw = zeros(720);
for m_min = 13:19
    m_max = m_min + 1;
    map = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
        num2str(m_min),'_',num2str(m_max),'sides.fits'));
    psmap_raw = map + psmap_raw;
end

%%%%%%%%%%%%%% mask %%%%%%%%%%%%%%
[mask,num]=make_mask_sides(flight,inst,0,20);

%%%%%%%%% mean/grad sub %%%%%%%%%%
grad_cb = plane_fit(cbmap_raw,mask);
grad_ps = plane_fit(psmap_raw,mask);
cb_bk = mean(cbmap_raw(find(mask)));
ps_bk = mean(psmap_raw(find(mask)));

cbmap = cbmap_raw - cb_bk - grad_cb;
psmap = psmap_raw - ps_bk - grad_ps;

stackmapdat(ifield).m_min_arr = [16:18];
stackmapdat(ifield).m_max_arr = [17:19];
stackmapdat(ifield).cbmap = cbmap;
stackmapdat(ifield).psmap = psmap;
stackmapdat(ifield).strmask = mask;
stackmapdat(ifield).strnum = num;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackmapdat_sides',savedir),'stackmapdat');

return