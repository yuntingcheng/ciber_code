function run_ihl_stack1(flight,inst,ifield)
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

%%% sig clip params %%%
nbins = 25;
iter_clip = 3;
sig = 2;
rmin = 20; % don't clip within rmin subpixels

for im=10%10:12%1:numel(m_min_arr)

%     if im==10
%         nbins = 28;
%         dx = 2000;
%     end
    
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    type = 1;
    [stampercb,stamperps,hitmap,count_totg,count_stackg]=...
        stackihl_ps(flight,inst,ifield,type,m_min,m_max,dx,...
        cbmap,psmap,mask_inst,strmask,strnum,nbins,rmin,iter_clip,sig,verbose);
    
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'_1.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmapg',...
        num2str(m_min),'_',num2str(m_max),'_1.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'_1.fits'),stamperps);

    type = -1;
    [stampercb,stamperps,hitmap,count_tots,count_stacks]=...
        stackihl_ps(flight,inst,ifield,type,m_min,m_max,dx,...
        cbmap,psmap,mask_inst,strmask,strnum,nbins,rmin,iter_clip,sig,verbose);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'_1.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'_1.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'_1.fits'),stamperps);
       
end

return