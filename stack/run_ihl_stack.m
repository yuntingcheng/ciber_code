function run_ihl_stack(flight,inst,ifield)
%%
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

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

if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end

for im=10:12
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    type = 1;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_ps_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum,stackband);

    [stampercbc,stamperpsc,hitmapc,idx_stack_arr]=...
        stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,0,verbose,0,0,stackband);
    
    [stampercb,stamperps,hitmap]=...
        stackihl_ps0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
        mask_inst,strmask,strnum,0,0,verbose,0,0,stackband,idx_stack_arr);
%     stampercb = fits_read(strcat(savedir,'ciber_ps/',dt.name,'_stamperg',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));
%     hitmap = load(strcat(savedir,'ciber_ps/',dt.name,'_hitmapg',...
%         num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
%     hitmap = hitmap.hitmap;
%     stamperps = fits_read(strcat(savedir,'panstarrs_ps/',dt.name,'_stamperg',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));

    stampercb(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stampercbc;
    stamperps(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stamperpsc;
    hitmap(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = hitmapc;
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmapg',...
        num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);

    stackcountdat(im).countg = numel(idx_stack_arr);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type = -1;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_ps_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum,stackband);

    [stampercbc,stamperpsc,hitmapc,idx_stack_arr]=...
        stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,0,verbose,0,0,stackband);
    
    [stampercb,stamperps,hitmap]=...
        stackihl_ps0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
        mask_inst,strmask,strnum,0,0,verbose,0,0,stackband,idx_stack_arr);
%     stampercb = fits_read(strcat(savedir,'ciber_ps/',dt.name,'_stampers',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));
%     hitmap = load(strcat(savedir,'ciber_ps/',dt.name,'_hitmaps',...
%         num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
%     hitmap = hitmap.hitmap;
%     stamperps = fits_read(strcat(savedir,'panstarrs_ps/',dt.name,'_stampers',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));

    stampercb(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stampercbc;
    stamperps(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stamperpsc;
    hitmap(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = hitmapc;
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
    
    stackcountdat(im).counts = numel(idx_stack_arr);
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackcountdat_%s',savedir,dt.name),'stackcountdat');

return