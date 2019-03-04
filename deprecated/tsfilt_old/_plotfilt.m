%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot the lab dark diff and flight diff PS filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% raw PS of flight and rail
flight=40030;
inst=1;
pixscale=7;
%ifield=5;
for ifield=[4,3,2,1]
savedir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
            'sinfiltamp/plot/field',num2str(ifield),'/');
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/DiffMap/field',num2str(ifield),'/');
load(strcat(loaddir,'maskin'),'maskin');
load(strcat(loaddir,'flightmap'),'flightmap');
dt=get_dark_times(flight,inst,ifield);
nfr_arr=2:dt.nfrhalf;
%%% get Fourier weight
rCl2dstack_arr=zeros(numel(nfr_arr),1024,1024);
fCl2dstack_arr=zeros(numel(nfr_arr),1024,1024);
for infr=1:numel(nfr_arr)
    infr
rCl2dstack=zeros(1024);
fCl2dstack=zeros(1024);    
for i=1:numel(dt.time)
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
end
rCl2dstack_arr(infr,:,:)=rCl2dstack;
fCl2dstack_arr(infr,:,:)=fCl2dstack;
end
save(strcat(loaddir,'fCl2dstack_arr'),'fCl2dstack_arr');
%%%
%%% plot Cl-l
rawdarkCl_arr=zeros(numel(dt.time),numel(nfr_arr),29);
filtdarkCl_arr=zeros(numel(dt.time),numel(nfr_arr),29);
rawflightCl_arr=zeros(1,numel(nfr_arr),29);
filtflightCl_arr=zeros(1,numel(nfr_arr),29);

for infr=1:numel(nfr_arr)  
    nfr=nfr_arr(infr);
    rCl2dstack=squeeze(rCl2dstack_arr(infr,:,:));
    fCl2dstack=squeeze(fCl2dstack_arr(infr,:,:));
    Wraw=(fftshift(fftshift(1./rCl2dstack)))';
    Wfilt=(fftshift(fftshift(1./fCl2dstack)))';
    figure 
    setwinsize(gcf,1000,500)
   %%%%%%%% raw %%%%%%%%
   subplot(1,2,1)
   for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rawmap=squeeze(labdat.rawmap_arr(infr,:,:));
        rawmap=rawmap-mean(rawmap(find(maskin)));
        rawmap=rawmap.*maskin;
        [labCl,l] = get_angular_spec(rawmap,rawmap,pixscale,'w',Wraw);
        rawdarkCl_arr(i,infr,:)=labCl;
        if i~=1
        loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end        
    end
    rawmap=squeeze(flightmap.rawmap_arr(infr,:,:));
    rawmap=rawmap-mean(rawmap(find(maskin)));
    rawmap=rawmap.*maskin;
    [rCl,l] = get_angular_spec(rawmap,rawmap,pixscale,'w',Wraw);
    rawflightCl_arr(1,infr,:)=rCl;
    pltflight=loglog(l,sqrt(l.*(l+1).*rCl),'ro','MarkerSize',5);hold off 
    xlim([1e2,2e5]);
    axlim=axis;
    legend([pltflight,pltdc],...
    {'flight unfilt','dark unfilt'},'location','southeast');
    legend boxoff
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    %%%%%%%% filt %%%%%%%%
    subplot(1,2,2)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        filtmap=squeeze(labdat.filtmap_arr(infr,:,:));
        filtmap=filtmap-mean(filtmap(find(maskin)));
        filtmap=filtmap.*maskin;
        [labCl,l] = get_angular_spec(filtmap,filtmap,pixscale,'w',Wfilt);
        filtdarkCl_arr(i,infr,:)=labCl;
        if i~=1
        loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=loglog(l,sqrt(l.*(l+1).*labCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end        
    end
    filtmap=squeeze(flightmap.filtmap_arr(infr,:,:));
    filtmap=filtmap-mean(filtmap(find(maskin)));
    filtmap=filtmap.*maskin;
    [fCl,l] = get_angular_spec(filtmap,filtmap,pixscale,'w',Wfilt);
    filtflightCl_arr(1,infr,:)=fCl;
    pltfunfilt=loglog(l,sqrt(l.*(l+1).*rCl),'r+','MarkerSize',5);hold on
    pltffilt=loglog(l,sqrt(l.*(l+1).*fCl),'ro','MarkerSize',5);hold off 
    xlim([axlim(1),axlim(2)]);
    ylim([axlim(3),axlim(4)]);
    legend([pltfunfilt,pltffilt,pltdc],...
    {'flight unfilt','filtght filt','dark filt'},'location','southeast');
    legend boxoff
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    savename=strcat(savedir,'nfr',num2str(nfr));
    print(savename,'-dpng');close
end
%%% plot Cl-nfr
for il=9:29
    lv=l(il);
    figure
    setwinsize(gcf,1000,500)
    
    %%%%%%%% raw %%%%%%%%
    subplot(1,2,1)
    for i=1:numel(dt.time)
        lCl=squeeze(rawdarkCl_arr(i,:,il));
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*lCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end  
    rawfCl=squeeze(rawflightCl_arr(1,:,il));
    pltf=semilogy(nfr_arr,sqrt(lv.*(lv+1).*rawfCl),'ro','MarkerSize',5);hold off
    axlim=axis;
    title(sprintf('unfilt,Fweight, l=%d',il))
    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    %%%%%%%% filt %%%%%%%%
    subplot(1,2,2)
    for i=1:numel(dt.time)
        lCl=squeeze(filtdarkCl_arr(i,:,il));
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*lCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end  
    filtfCl=squeeze(filtflightCl_arr(1,:,il));
    pltf=semilogy(nfr_arr,sqrt(lv.*(lv+1).*filtfCl),'ro','MarkerSize',5);hold off
    xlim([axlim(1),axlim(2)]);
    ylim([axlim(3),axlim(4)]);
    title(sprintf('filt,Fweight, l=%d',il))
    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    savename=strcat(savedir,'ell',num2str(il));
    print(savename,'-dpng');close    
end
%%%
save(strcat(savedir,'rawdarkCl_arr'),'rawdarkCl_arr');
save(strcat(savedir,'rawflightCl_arr'),'rawflightCl_arr');
save(strcat(savedir,'filtdarkCl_arr'),'filtdarkCl_arr');
save(strcat(savedir,'filtflightCl_arr'),'filtflightCl_arr');
end