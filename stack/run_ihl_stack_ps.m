flight=40030;
inst=1;
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));
m_min_arr = [0,8:22];
m_max_arr = [8:23];
%% get CIBER and 2MASS map

for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    %%%%% CIBER map (from make_lin_premap.m) %%%%%%%
    cbmap_raw = fits_read(strcat(savedir,'maps/',dt.name,'_map.fits'));
    % correct for the cal factor error
    cbmap_raw = cbmap_raw./0.383;

    %%%%% PanSTARRS srcmap %%%%%%%
    srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
        num2str(inst),'/');
    psmap_raw = zeros(1024);
    for im=1:numel(m_min_arr)
        psmapsi = fits_read(strcat(srcmapdir,dt.name,'_srcmaps',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps.fits'));
        psmapgi = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps.fits'));
        psmapui = fits_read(strcat(srcmapdir,dt.name,'_srcmapu',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps.fits'));
        psmap_raw = psmap_raw + psmapsi + psmapgi + psmapui;
    end

    %%% masks %%%
    mask_inst = fits_read(strcat(savedir,'masks/',dt.name,'_mask.fits'));
    strmask = fits_read(strcat(savedir,'masksps/',dt.name,'_strmask_all.fits'));
    strnum = fits_read(strcat(savedir,'masksps/',dt.name,'_strnum_all.fits'));
    totmask = mask_inst.*strmask;

    %%% sigma clip and mean/grad sub %%%
    sigmask = sigclip_mask(cbmap_raw,totmask,5,5);
    grad_cb = plane_fit(cbmap_raw,sigmask);
    grad_ps = plane_fit(psmap_raw,sigmask);
    cb_bk = mean(cbmap_raw(find(sigmask)));
    ps_bk = mean(psmap_raw(find(sigmask)));

    cbmap = cbmap_raw - cb_bk - grad_cb;
    psmap = psmap_raw - ps_bk - grad_ps;

    sig_sp = find((totmask-sigmask)==1);
    mask_inst_clip = mask_inst;
    mask_inst_clip(sig_sp)=0;

    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_cbmap.fits'),...
        cbmap);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_psmap.fits'),...
        psmap);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_mask_inst_clip.fits'),...
        mask_inst_clip);
end

%% stacking function

dx = 1200;
verbose=1;
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    
    type = -1;
    [stampercb,stamperps,maskstamper,count_tots,count_stacks]=...
        stackihl_ps(flight,inst,ifield,type,m_min,m_max,dx,...
        cbmap,psmap,mask_inst_clip,strmask,strnum,verbose);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'.fits'),maskstamper);
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
    
    
    type = 1;
    [stampercb,stamperps,maskstamper,count_totg,count_stackg]=...
        stackihl_ps(flight,inst,ifield,type,m_min,m_max,dx,...
        cbmap,psmap,mask_inst_clip,strmask,strnum,verbose);
    
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_hitmapg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),maskstamper);
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
        
end
end