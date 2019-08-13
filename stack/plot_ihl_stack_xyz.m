flight=40030;
inst=1;
mypaths=get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

for ifield=4:8
spire=false;
rmin=nan;
mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if rmin==2
    load(sprintf('%sstackmapdat_rmin2',loaddir),'stackmapdat');
elseif isnan(rmin)
    load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
end
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
return_count=true;

if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end
%%
figure
setwinsize(gcf,1200,750)

for im= 10:13
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    type = -1;
    [~,~,~,~,stackcount,subx_arr,suby_arr,~]=...
    stackihl_ps0_hist(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,1,0,verbose,10,stackband,[],spire,rmin,return_count);    
    hitmaps = zeros(1024);
    for i=1:stackcount
        hitmaps(round(subx_arr(i)),round(suby_arr(i))) = ...
            hitmaps(round(subx_arr(i)),round(suby_arr(i))) + 1;
    end
       
    type = 1;
    [~,~,~,~,stackcount,subx_arr,suby_arr,subz_arr]=...
    stackihl_ps0_hist(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,1,0,verbose,10,stackband,[],spire,rmin,return_count);
    hitmapg = zeros(1024);
    for i=1:stackcount
        hitmapg(round(subx_arr(i)),round(suby_arr(i))) = ...
            hitmapg(round(subx_arr(i)),round(suby_arr(i))) + 1;
    end
    
    subplot(3,4,im-9)
    imageclip(rebin_map_coarse(hitmaps,8).*64);
    caxis([0,3]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    title(strcat(num2str(m_min),'<m<',num2str(m_max),' stars'),'fontsize',15);

    subplot(3,4,im-9+4)
    imageclip(rebin_map_coarse(hitmapg,8).*64);
    caxis([0,3]);
    set(gca,'XTickLabel',[]);
    set(gca,'YTickLabel',[]);
    title(strcat(num2str(m_min),'<m<',num2str(m_max),' galaxies'),'fontsize',15);
    
    subplot(3,4,im-9+8)
    binedges = 0:0.05:1.1;
    h = histcounts(subz_arr,'BinEdges',binedges) / (binedges(2)-binedges(1));
    bins = binedges2bins(binedges);
    plot(bins,h,'LineWidth',3);
    xlim([0,1.1]);
    text(sum(get(gca,'XLim').*[0.6,0.4]), sum(get(gca,'YLim').*[0.1,0.9]),...
    sprintf('median z = %.2f', median(subz_arr)),'fontsize',15);
    xlabel('z','fontsize',15);
    ylabel('dN/dz','fontsize',15);
    title('stacked galaxy redshift','fontsize',12);
end
s=suptitle(dt.name);
set(s,'FontSize',15);
savename=strcat(pltsavedir,dt.name,'_stackxyz');
print(savename,'-dpng');
end

