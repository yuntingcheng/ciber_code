function plot_filt_2DPS(flight,inst,ifield)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the individual RN 2DCl, stacked RN 2DCl (1st-2nd diff)
%before and after filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pixscale=7;
mypaths=get_paths(flight);

loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
load(strcat(loaddir,'bigmask'),'bigmask');
dt=get_dark_times(flight,inst,ifield);

ell=get_l(1024,1024,pixscale,1);
rCl2dstack=zeros(1024);
fCl2dstack=zeros(1024);
for i=1:numel(dt.time)
    load(strcat(loaddir,'labmap',num2str(i)),'labmap');
    rawmap=(labmap.rawmap1-labmap.rawmap2)./2;
    [~,maskin1]=get_skymap(rawmap,bigmask,4,5);
    rawmap=rawmap-mean(rawmap(find(maskin1)));
    rawmap=rawmap.*maskin1;
    [rCl,l,~,~,binl,~,rCl2d] = get_angular_spec(rawmap,rawmap,pixscale);
    rCl2dstack=rCl2dstack+rCl2d./numel(dt.time);
    
    filtmap=(labmap.filtmap1-labmap.filtmap2)./2;
    [~,maskin1]=get_skymap(filtmap,bigmask,4,5);
    filtmap=filtmap-mean(filtmap(find(maskin1)));
    filtmap=filtmap.*maskin1;
    [fCl,l,~,~,~,~,fCl2d] = get_angular_spec(filtmap,filtmap,pixscale);
    fCl2dstack=fCl2dstack+fCl2d./numel(dt.time);   
    
    fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
    [x,y]=find(fmask);
    
    figure
    setwinsize(gcf,1000,800)
    
    subplot(2,2,1)
    imageclip(rawmap);
    v=caxis;
    title('unfilt map');
    
    subplot(2,2,2)
    imageclip(filtmap);
    caxis(v);
    title('filt map');
    
    subplot(2,2,3)
    imageclip(rCl2d);
    v=caxis;
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    title('unfilt 2D Cl');
    
    subplot(2,2,4)
    imageclip(fCl2d);
    caxis(v);
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    title('filt 2D Cl');
    
    drawnow
    %savename=strcat(plotdir,'mapCl_single');
    %print(savename,'-dpng');close    
end

%%%%%%% now plot thestacked 2DPS %%%%%%%%%%
fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
[x,y]=find(fmask);

figure
setwinsize(gcf,1000,800)
subplot(2,2,1)
imageclip(rCl2dstack);
v=caxis;
title('unfilt');

subplot(2,2,2)
imageclip(fCl2dstack);
caxis(v);
title('filt');

subplot(2,2,3)
imageclip(rCl2dstack(min(x):max(x),min(y):max(y)));
v=caxis;

subplot(2,2,4)
imageclip(fCl2dstack(min(x):max(x),min(y):max(y)));
caxis(v);

%savename=strcat(plotdir,'mapCl_stack');
%print(savename,'-dpng');close
end