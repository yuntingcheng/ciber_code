function run_ihl_stack2(flight,inst,ifield)
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

nbins = 25;
rmin = 15;
sig = 20;
%%
for im=12%10:12
%%
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    maskrad = get_mask_radius(inst,ifield,(m_min + m_max)/2);
    profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
    r_arr=profile.r*0.7;

    cliplimg_arr = ihl_stack_cliplim(flight,inst,ifield,1,...
    m_max,m_min,dx,nbins,sig,rmin);
    cliplims_arr = ihl_stack_cliplim(flight,inst,ifield,-1,...
    m_max,m_min,dx,nbins,sig,rmin);
    
    figure
    loglog(r_arr,cliplimg_arr(1,:),'bo-');hold on
    loglog(r_arr,cliplims_arr(1,:),'ro-');
    loglog(r_arr,cliplimg_arr(3,:),'co-');
    loglog(r_arr,cliplims_arr(3,:),'mo-');

    for i=1:4
        means = median(cliplims_arr(i,r_arr > maskrad));
        meang = median(cliplimg_arr(i,r_arr > maskrad));
        if im==12
            cliplimg_arr(i,r_arr > maskrad) = min([means,meang]);
            cliplims_arr(i,r_arr > maskrad) = min([means,meang]);
        else
            cliplimg_arr(i,r_arr > maskrad) = meang;
            cliplims_arr(i,r_arr > maskrad) = means;            
        end
    end
    loglog(r_arr,cliplimg_arr(1,:),'b.-');hold on
    loglog(r_arr,cliplims_arr(1,:),'r.-');
    loglog(r_arr,cliplimg_arr(3,:),'c.-');
    loglog(r_arr,cliplims_arr(3,:),'m.-');
    
    type = 1;
    [stampercb,stamperps,hitmap,count_totg,count_stackg]=...
        stackihl_ps2(flight,inst,ifield,type,m_min,m_max,dx,...
        cbmap,psmap,mask_inst,strmask,strnum,nbins,rmin,cliplimg_arr,verbose);
%%
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmapg',...
        num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stamperg',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);

    type = -1;
    [stampercb,stamperps,hitmap,count_tots,count_stacks]=...
        stackihl_ps2(flight,inst,ifield,type,m_min,m_max,dx,...
        cbmap,psmap,mask_inst,strmask,strnum,nbins,rmin,cliplims_arr,verbose);
    fits_write(strcat(savedir,'ciber_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stampercb);
    save(strcat(savedir,'ciber_ps/',dt.name,'_hitmaps',...
        num2str(m_min),'_',num2str(m_max),'.mat'),'hitmap');
    fits_write(strcat(savedir,'panstarrs_ps/',dt.name,'_stampers',...
        num2str(m_min),'_',num2str(m_max),'.fits'),stamperps);
end

return