function plot_Cl_ph(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot filtered diff RN + ph noise w/ G1 fitted in
% individual and all fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

pixscale=7;
cp=get_cal_params('flight',flight);
frate=cp(inst).framerate;

%%% load best g1
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/');
load(strcat(loaddir,'diffG1fit/chi2dat'),'chi2dat');

%%% load fligt spectrum
load(sprintf('%s/diffCldat',loaddir),'diffCldat');

for ifield=1:8
G1tot=chi2dat.bestg1_allfields;
G1=chi2dat.bestg1(ifield).g1;
if isempty(G1);G1=G1tot;end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),...
    '/noisemodel/diffplotph/field',num2str(ifield),'/');
dt=get_dark_times(flight,inst,ifield);
loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
load(strcat(loaddir,'bigmask'),'bigmask');
load(strcat(loaddir,'flightmap'),'flightmap');

nfr_arr=2:dt.nfrhalf;
slopemap=flightmap.rawmapf;
slopemean=mean(slopemap(find(bigmask)));

%%% plot Cl-l
nClsingG1_arr=zeros(numel(dt.time),numel(nfr_arr),29);
nCltotG1_arr=zeros(numel(dt.time),numel(nfr_arr),29);

for infr=1:numel(nfr_arr)  
    nfr=nfr_arr(infr);
    %%% get fweight %%%
    fCl2d_arr=zeros(numel(dt.time),1024,1024);
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');
        filtmap=squeeze(labmap.filtmap_arr(nfr-1,:,:));
        [~,maskin1]=get_skymap(filtmap,bigmask,4,5);
        filtmap=filtmap-mean(filtmap(find(maskin1)));
        filtmap=filtmap.*maskin1;
        filtmap=dc_offset_remove(filtmap,maskin1).*maskin1;
        [~,~,~,~,~,~,fCl2d] = get_angular_spec(filtmap,filtmap,pixscale);
        fCl2d_arr(i,:,:)=fCl2d;
    end
    fCl2d_std=squeeze(std(fCl2d_arr));
    fw=(fftshift(fftshift(1./squeeze(fCl2d_std))))';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure
    setwinsize(gcf,1000,500)
    %%%%%%%% indiv G1 %%%%%%%%
    subplot(1,2,1)
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');
        rnmap=squeeze(labmap.filtmap_arr(infr,:,:));        
        phmap1=photonnoise_realization(ones(1024).*slopemean,G1,nfr,frate);
        phmap2=photonnoise_realization(ones(1024).*slopemean,G1,nfr,frate);
        phmap=(phmap1-phmap2)./2;

        nmap=rnmap+phmap;
        [~,maskin1]=get_skymap(nmap,bigmask,4,5);
        nmap=nmap-mean(nmap(find(maskin1)));
        nmap=nmap.*maskin1;
        [nCl,l] = get_angular_spec(nmap,nmap,pixscale,'w',fw);
        nClsingG1_arr(i,infr,:)=nCl;
        if i~=1
        loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end        
    end
    rCl=squeeze(diffCldat.flight(ifield).wrCl_arr(infr,:));
    fCl=squeeze(diffCldat.flight(ifield).wfCl_arr(infr,:));
    pltfunfilt=loglog(l,sqrt(l.*(l+1).*rCl),'r+','MarkerSize',5);hold on
    pltffilt=loglog(l,sqrt(l.*(l+1).*fCl),'ro','MarkerSize',5);hold off 
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
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');
        rnmap=squeeze(labmap.filtmap_arr(infr,:,:));        
        phmap1=photonnoise_realization(ones(1024).*slopemean,G1tot,nfr,frate);
        phmap2=photonnoise_realization(ones(1024).*slopemean,G1tot,nfr,frate);
        phmap=(phmap1-phmap2)./2;

        nmap=rnmap+phmap;
        [~,maskin1]=get_skymap(nmap,bigmask,4,5);
        nmap=nmap-mean(nmap(find(maskin1)));
        nmap=nmap.*maskin1;
        [nCl,l] = get_angular_spec(nmap,nmap,pixscale,'w',fw);
        nCltotG1_arr(i,infr,:)=nCl;
        if i~=1
        loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        else
        pltdc=loglog(l,sqrt(l.*(l+1).*nCl),'color',[0.7,0.7,0.7]);hold on
        drawnow
        end        
    end
    rCl=squeeze(diffCldat.flight(ifield).wrCl_arr(infr,:));
    fCl=squeeze(diffCldat.flight(ifield).wfCl_arr(infr,:));
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
    title(sprintf('best fit all fields, G1=%.2f',G1tot));
    drawnow
    
    savename=strcat(savedir,'nfr',num2str(nfr));
    print(savename,'-dpng');close
end
diffCldat.dark(ifield).nG1indivCl_arr=nClsingG1_arr;
diffCldat.dark(ifield).nG1Cl_arr=nCltotG1_arr;
end
save(sprintf('%sTM%d/diffCldat',mypaths.filtmap,inst),'diffCldat');

%%%%%%%%%%%%%% plot Cl-nfr %%%%%%%%%%%%%%%%%%

load(sprintf('%sTM%d/diffCldat',mypaths.filtmap,inst),'diffCldat');
for ifield=4:8%%%
savedir=strcat(mypaths.alldat,'TM',num2str(inst),...
    '/noisemodel/diffplotph/field',num2str(ifield),'/');
dt=get_dark_times(flight,inst,ifield);
nfr_arr=2:dt.nfrhalf;
G1tot=chi2dat.bestg1_allfields;
G1=chi2dat.bestg1(ifield).g1;
nClsingG1_arr=diffCldat.dark(ifield).nG1indivCl_arr;
nCltotG1_arr=diffCldat.dark(ifield).nG1Cl_arr;
[~,l]=get_angular_spec(randn(1024),randn(1024),pixscale);
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
    rCl=squeeze(diffCldat.flight(ifield).wrCl_arr(:,il));
    fCl=squeeze(diffCldat.flight(ifield).wfCl_arr(:,il));
    pltr=semilogy(nfr_arr,sqrt(lv.*(lv+1).*rCl),'r+','MarkerSize',5);hold on
    pltf=semilogy(nfr_arr,sqrt(lv.*(lv+1).*fCl),'ro','MarkerSize',5);hold off
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
    rCl=squeeze(diffCldat.flight(ifield).wrCl_arr(:,il));
    fCl=squeeze(diffCldat.flight(ifield).wfCl_arr(:,il));
    pltr=semilogy(nfr_arr,sqrt(lv.*(lv+1).*rCl),'r+','MarkerSize',5);hold on
    pltf=semilogy(nfr_arr,sqrt(lv.*(lv+1).*fCl),'ro','MarkerSize',5);hold off
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

return