%% raw PS of flight and rail
flight=40030;
inst=1;
pixscale=7;
ifield=8;

savedir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
            'sinfiltrange/plot/field',num2str(ifield),'/');
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltrange/DiffMap/field',num2str(ifield),'/');
plotdir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
       'sinfiltrange/plot/field',num2str(ifield),'/');

load(strcat(loaddir,'maskin'),'maskin');
load(strcat(loaddir,'flightmap'),'flightmap');
dt=get_dark_times(flight,inst,ifield);
nfr_arr=2:dt.nfrhalf;
%%
ell=get_l(1024,1024,pixscale,1);
rCl2dstack=zeros(1024);
fCl2dstack=zeros(1024);
infr=12;
for i=1:numel(dt.time)
    i
    load(strcat(loaddir,'labdat',num2str(i)),'labdat');
    rawmap=squeeze(labdat.rawmap_arr(infr,:,:));
    rawmap=rawmap-mean(rawmap(find(maskin)));
    rawmap=rawmap.*maskin;
    [rCl,l,~,~,binl,~,rCl2d] = get_angular_spec(rawmap,rawmap,pixscale);
    rCl2dstack=rCl2dstack+rCl2d./numel(dt.time);
    
    filtmap=squeeze(labdat.filtmap_arr(infr,:,:));
    filtmap=filtmap-mean(filtmap(find(maskin)));
    filtmap=filtmap.*maskin;
    [fCl,l,~,~,~,~,fCl2d] = get_angular_spec(filtmap,filtmap,pixscale);
    fCl2dstack=fCl2dstack+fCl2d./numel(dt.time);   
    
    fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
    [x,y]=find(fmask);
    
    if i==1
    figure
    setwinsize(gcf,1000,800)
    subplot(2,2,1)
    imageclip(rCl2d);
    v=caxis;
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    title('unfilt 2D Cl');
    
    subplot(2,2,2)
    imageclip(fCl2d);
    caxis(v);
    xlim([min(x),max(x)]); ylim([min(y),max(y)]);
    title('filt 2D Cl');
    
    subplot(2,2,3)
    imageclip(rawmap);
    v=caxis;
    title('unfilt map');
    
    subplot(2,2,4)
    imageclip(filtmap);
    caxis(v);
    title('filt map');
    
    drawnow
    savename=strcat(plotdir,'mapCl_single');
    print(savename,'-dpng');close

    end
    
end
%%
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

savename=strcat(plotdir,'mapCl_stack');
print(savename,'-dpng');close
%% plot unfilt
Wraw=(fftshift(fftshift(1./rCl2dstack)))';
for infr=12%1:numel(nfr_arr)  
    nfr=nfr_arr(infr);
    figure 
    setwinsize(gcf,1000,400)
        %%%%%%%% weight %%%%%%%%
    subplot(1,2,1)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rawmap=squeeze(labdat.rawmap_arr(infr,:,:));
        rawmap=rawmap-mean(rawmap(find(maskin)));
        rawmap=rawmap.*maskin;
        [labCl,l] = get_angular_spec(rawmap,rawmap,pixscale);
        loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end
    rawmap=squeeze(flightmap.rawmap_arr(infr,:,:));
    rawmap=rawmap-mean(rawmap(find(maskin)));
    rawmap=rawmap.*maskin;
    [fCl] = get_angular_spec(rawmap,rawmap,pixscale);

    pltf=loglog(l,sqrt(l.*(l+1).*fCl),'ro','MarkerSize',5);%hold off    
    title(sprintf('N=%d, unfiltered; unweighted',nfr))
    xlim([1e2,2e5]);
    axlim=axis;
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    %%%%%%%% weight %%%%%%%%
    subplot(1,2,2)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rawmap=squeeze(labdat.rawmap_arr(infr,:,:));
        rawmap=rawmap-mean(rawmap(find(maskin)));
        rawmap=rawmap.*maskin;
        [labCl,l] = get_angular_spec(rawmap,rawmap,pixscale,'w',Wraw);
        loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end
    rawmap=squeeze(flightmap.rawmap_arr(infr,:,:));
    rawmap=rawmap-mean(rawmap(find(maskin)));
    rawmap=rawmap.*maskin;
    [fCl] = get_angular_spec(rawmap,rawmap,pixscale,'w',Wraw);

    pltf=loglog(l,sqrt(l.*(l+1).*fCl),'ro','MarkerSize',5);%hold off    
    title(sprintf('N=%d, unfiltered; weighted',nfr))
    xlim([1e2,2e5]);
    axlim=axis;
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    savename=strcat(savedir,'unfilt_nfr',num2str(nfr));
    print(savename,'-dpng');close

end

%% plot filt
Wfilt=(fftshift(fftshift(1./fCl2dstack)))';
for infr=12%1:numel(nfr_arr)  
    nfr=nfr_arr(infr);
    figure 
    setwinsize(gcf,1000,400)
        %%%%%%%% weight %%%%%%%%
    subplot(1,2,1)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rawmap=squeeze(labdat.filtmap_arr(infr,:,:));
        rawmap=rawmap-mean(rawmap(find(maskin)));
        rawmap=rawmap.*maskin;
        [labCl,l] = get_angular_spec(rawmap,rawmap,pixscale);
        loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end
    rawmap=squeeze(flightmap.filtmap_arr(infr,:,:));
    rawmap=rawmap-mean(rawmap(find(maskin)));
    rawmap=rawmap.*maskin;
    [fCl] = get_angular_spec(rawmap,rawmap,pixscale);

    pltf=loglog(l,sqrt(l.*(l+1).*fCl),'ro','MarkerSize',5);%hold off    
    title(sprintf('N=%d, filtered; unweighted',nfr))
    xlim([1e2,2e5]);
    axlim=axis;
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    %%%%%%%% weight %%%%%%%%
    subplot(1,2,2)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rawmap=squeeze(labdat.filtmap_arr(infr,:,:));
        rawmap=rawmap-mean(filtmap(find(maskin)));
        rawmap=rawmap.*maskin;
        [labCl,l] = get_angular_spec(rawmap,rawmap,pixscale,'w',Wfilt);
        loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end
    rawmap=squeeze(flightmap.filtmap_arr(infr,:,:));
    rawmap=rawmap-mean(rawmap(find(maskin)));
    rawmap=rawmap.*maskin;
    [fCl] = get_angular_spec(rawmap,rawmap,pixscale,'w',Wfilt);

    pltf=loglog(l,sqrt(l.*(l+1).*fCl),'ro','MarkerSize',5);%hold off    
    title(sprintf('N=%d, filtered; weighted',nfr))
    xlim([1e2,2e5]);
    axlim=axis;
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    savename=strcat(savedir,'filt_nfr',num2str(nfr));
    print(savename,'-dpng');close

end