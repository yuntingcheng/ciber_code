%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the diff PS with best fit Chi2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% raw PS of flight and rail
flight=40030;
inst=1;
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;

%%% load best g1
bestg1dir='/Users/ytcheng/ciber/doc/20170209_TsFilter/G1fit/';
%load(strcat(bestg1dir,'bestg1'),'bestg1');
load(strcat(bestg1dir,'bestg1w'),'bestg1');
%%
for ifield=[8,7,6,5,4,2,1]
G1=bestg1(ifield);
G1tot=bestg1(9);
%savedir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
%            'G1fit/field',num2str(ifield),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
            'G1fit/whiteph/field',num2str(ifield),'/');

loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/DiffMap/field',num2str(ifield),'/');
load(strcat(loaddir,'maskin'),'maskin');
load(strcat(loaddir,'flightmap'),'flightmap');
dt=get_dark_times(flight,inst,ifield);
nfr_arr=2:dt.nfrhalf;
%%% get Fourier weight
load(strcat(loaddir,'fCl2dstack_arr'),'fCl2dstack_arr');
%%% load flight PS
loaddir1=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
            'sinfiltamp/plot/field',num2str(ifield),'/');
load(strcat(loaddir1,'rawflightCl_arr'),'rawflightCl_arr');
load(strcat(loaddir1,'filtflightCl_arr'),'filtflightCl_arr');

%%% get slope map
[flightfr] = get_data_frames...
        (inst,dt.name,'flight',flight,'verbose',0);
flightfr=flightfr(3:end,:,:);
%%%%%%%%
slopemap=linfit_map(flightfr,'verbose',0);
slopemean=mean(slopemap(find(maskin)));
slopemap1=ones(1024)*slopemean;
slopemap2=ones(1024)*slopemean;
%%%%%%%%
%%
%%% plot Cl-l
nClsingG1_arr=zeros(numel(dt.time),numel(nfr_arr),29);
nCltotG1_arr=zeros(numel(dt.time),numel(nfr_arr),29);

for infr=1:numel(nfr_arr)  
    nfr=nfr_arr(infr);
    fCl2dstack=squeeze(fCl2dstack_arr(infr,:,:));
    Wfilt=(fftshift(fftshift(1./fCl2dstack)))';
    %%%%%%%%%%%
    %fr1=flightfr(1:nfr,:,:);
    %fr2=flightfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);
    %[slopemap1]=linfit_map(fr1,'verbose',0);
    %[slopemap2]=linfit_map(fr2,'verbose',0);
    %%%%%%%%%%%%%%%
    figure
    setwinsize(gcf,1000,500)
    %%%%%%%% indiv G1 %%%%%%%%
    subplot(1,2,1)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rnmap=squeeze(labdat.filtmap_arr(infr,:,:));

        phmap1=photonnoise_realization(slopemap1,G1,nfr,frate);
        phmap2=photonnoise_realization(slopemap2,G1,nfr,frate);
        phmap=(phmap1-phmap2)./2;

        nmap=rnmap+phmap;
        nmap=nmap-mean(nmap(find(maskin)));
        nmap=nmap.*maskin;
        [nCl,l] = get_angular_spec(nmap,nmap,pixscale,'w',Wfilt);
        nClsingG1_arr(i,infr,:)=nCl;
        if i~=1
        loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end        
    end
    fCl=squeeze(filtflightCl_arr(1,infr,:));
    rCl=squeeze(rawflightCl_arr(1,infr,:));
    pltfunfilt=loglog(l,sqrt(l.*(l+1).*rCl'),'r+','MarkerSize',5);hold on
    pltffilt=loglog(l,sqrt(l.*(l+1).*fCl'),'ro','MarkerSize',5);hold off 
    xlim([1e2,2e5]);
    axlim=axis;
    legend([pltfunfilt,pltffilt,pltdc],...
    {'flight unfilt','filtght filt','dark filt+photon noise'},...
    'location','southeast');
    legend boxoff
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    title(sprintf('best fit single field, G1=%.2f',G1));
    drawnow
    %%%%%%%% tot G1 %%%%%%%%
    subplot(1,2,2)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rnmap=squeeze(labdat.filtmap_arr(infr,:,:));

        phmap1=photonnoise_realization(slopemap1,G1tot,nfr,frate);
        phmap2=photonnoise_realization(slopemap2,G1tot,nfr,frate);
        phmap=(phmap1-phmap2)./2;

        nmap=rnmap+phmap;
        nmap=nmap-mean(nmap(find(maskin)));
        nmap=nmap.*maskin;
        [nCl,l] = get_angular_spec(nmap,nmap,pixscale,'w',Wfilt);
        nCltotG1_arr(i,infr,:)=nCl;
        if i~=1
        loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end        
    end
    fCl=squeeze(filtflightCl_arr(1,infr,:));
    rCl=squeeze(rawflightCl_arr(1,infr,:));
    pltfunfilt=loglog(l,sqrt(l.*(l+1).*rCl'),'r+','MarkerSize',5);hold on
    pltffilt=loglog(l,sqrt(l.*(l+1).*fCl'),'ro','MarkerSize',5);hold off 
    xlim([axlim(1),axlim(2)]);
    ylim([axlim(3),axlim(4)]);
    legend([pltfunfilt,pltffilt,pltdc],...
    {'flight unfilt','filtght filt','dark filt'},'location','southeast');
    legend boxoff
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    title(sprintf('best fit all fields, G1=%.2f',G1tot));
    drawnow
    
    savename=strcat(savedir,'nfr',num2str(nfr));
    print(savename,'-dpng');close
end
save(strcat(savedir,'nClsingG1_arr'),'nClsingG1_arr');
save(strcat(savedir,'nCltotG1_arr'),'nCltotG1_arr');




%%% plot Cl-nfr
for il=9:29
    lv=l(il);
    figure
    setwinsize(gcf,1000,500)
    
    %%%%%%%% indiv G1 %%%%%%%%
    subplot(1,2,1)
    for i=1:numel(dt.time)
        lCl=squeeze(nClsingG1_arr(i,:,il));
        if i~=1
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*lCl),'color',[0.7,0.7,0.7]);hold on
        else
        pltdc=semilogy(nfr_arr,sqrt(lv.*(lv+1).*lCl),'color',[0.7,0.7,0.7]);hold on
        end
       
        drawnow
    end  
    rawfCl=squeeze(rawflightCl_arr(1,:,il));
    pltr=semilogy(nfr_arr,sqrt(lv.*(lv+1).*rawfCl),'r+','MarkerSize',5);hold on
    filtfCl=squeeze(filtflightCl_arr(1,:,il));
    pltf=semilogy(nfr_arr,sqrt(lv.*(lv+1).*filtfCl),'ro','MarkerSize',5);hold off
    axlim=axis;
    title(sprintf('best fit single field, G1=%.2f',G1));
    legend([pltr,pltf,pltdc],...
    {'flight unfilt','filtght filt','dark filt'},'location','northeast');
    legend boxoff

    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    %%%%%%%% tot G1 %%%%%%%%
    subplot(1,2,2)
    for i=1:numel(dt.time)
        lCl=squeeze(nCltotG1_arr(i,:,il));
        if i~=1
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*lCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=semilogy(nfr_arr,sqrt(lv.*(lv+1).*lCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end

    end 
    rawfCl=squeeze(rawflightCl_arr(1,:,il));
    pltr=semilogy(nfr_arr,sqrt(lv.*(lv+1).*rawfCl),'r+','MarkerSize',5);hold on
    filtfCl=squeeze(filtflightCl_arr(1,:,il));
    pltf=semilogy(nfr_arr,sqrt(lv.*(lv+1).*filtfCl),'ro','MarkerSize',5);hold off
    xlim([axlim(1),axlim(2)]);
    ylim([axlim(3),axlim(4)]);
    title(sprintf('best fit all fields, G1=%.2f',G1tot));
    legend([pltr,pltf,pltdc],...
    {'flight unfilt','filtght filt','dark filt'},'location','northeast');
    legend boxoff

    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    savename=strcat(savedir,'ell',num2str(il));
    print(savename,'-dpng');close    
end
end