function run_ihl_stack_sim(flight,inst,ifield,set)
% set: 0 - all sim sources in CBmap, 
%      1 - m<22 sources in CBmap
%%
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdatsim',loaddir),'stackmapdatsim');

dx = 1200;
verbose = 1;
dt=get_dark_times(flight,inst,ifield);
if set==0
    cbmap = stackmapdatsim(ifield).all.cbmap;
    psmap = stackmapdatsim(ifield).all.psmap;
    mask_inst = stackmapdatsim(ifield).all.mask_inst_clip;
    strmask = stackmapdatsim(ifield).all.strmask;
    strnum = stackmapdatsim(ifield).all.strnum;
elseif set==1
    cbmap = stackmapdatsim(ifield).sub.cbmap;
    psmap = stackmapdatsim(ifield).sub.psmap;
    mask_inst = stackmapdatsim(ifield).sub.mask_inst_clip;
    strmask = stackmapdatsim(ifield).sub.strmask;
    strnum = stackmapdatsim(ifield).sub.strnum;
end

for im=1:3
    m_min = im + 15;
    m_max = m_min + 1;

    type = 1;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_sim_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum);

    [stampercbc,stamperpsc,hitmapc,idx_stack_arr]=...
        stackihl_sim_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,0,verbose,0,0);
    
    [stampercb,stamperps,hitmap]=...
        stackihl_sim0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
        mask_inst,strmask,strnum,0,0,verbose,0,0,idx_stack_arr);
%     stampercb = fits_read(strcat(savedir,'ciber_sim/',dt.name,'_stamperg',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));
%     hitmap = load(strcat(savedir,'ciber_sim/',dt.name,'_hitmapg',...
%         num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
%     hitmap = hitmap.hitmap;
%     stamperps = fits_read(strcat(savedir,'panstarrs_sim/',dt.name,'_stamperg',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));
    stampercb(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stampercbc;
    stamperps(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stamperpsc;
    hitmap(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = hitmapc;
    fits_write(strcat(savedir,'ciber_sim/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'),stampercb);
    save(strcat(savedir,'ciber_sim/',dt.name,'_hitmapg',...
        num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_sim/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'),stamperps);

    stackcountdat(im).countg = numel(idx_stack_arr);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    type = -1;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_sim_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum);

    [stampercbc,stamperpsc,hitmapc,idx_stack_arr]=...
        stackihl_sim_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,0,verbose,0,0);
    
    [stampercb,stamperps,hitmap]=...
        stackihl_sim0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
        mask_inst,strmask,strnum,0,0,verbose,0,0,idx_stack_arr);
%     stampercb = fits_read(strcat(savedir,'ciber_sim/',dt.name,'_stampers',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));
%     hitmap = load(strcat(savedir,'ciber_sim/',dt.name,'_hitmaps',...
%         num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
%     hitmap = hitmap.hitmap;
%     stamperps = fits_read(strcat(savedir,'panstarrs_sim/',dt.name,'_stampers',...
%         num2str(m_min),'_',num2str(m_max),'.fits'));
    stampercb(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stampercbc;
    stamperps(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stamperpsc;
    hitmap(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = hitmapc;
    fits_write(strcat(savedir,'ciber_sim/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'),stampercb);
    save(strcat(savedir,'ciber_sim/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_sim/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'_',num2str(set),'.fits'),stamperps);
    
    stackcountdat(im).counts = numel(idx_stack_arr);
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackcountdatsim_%s_%d',savedir,dt.name,set),'stackcountdat');

return