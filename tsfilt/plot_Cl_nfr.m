function plot_Cl_nfr(flight,inst)
mypaths=get_paths(flight);

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/diffplot/');
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/');
load(sprintf('%s/diffCldat',loaddir),'diffCldat');

l=diffCldat.l;
for ifield=1:8
dt=get_dark_times(flight,inst,ifield);
nfr_arr=diffCldat.dark(ifield).nfr_arr;
for il=9:29
    lv=l(il);
    
    figure
    setwinsize(gcf,1000,1000)
    
    subplot(2,2,1)
    for i=1:numel(dt.time)
        Cl=squeeze(diffCldat.dark(ifield).rCl_arr(i,:,il));
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end  
    Cl=squeeze(diffCldat.flight(ifield).rCl_arr(:,il))';
    semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'r+','MarkerSize',3);hold off
    axlim=axis;
    title(sprintf('%s, l=%d',dt.name,il))
    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    
    subplot(2,2,2)
    for i=1:numel(dt.time)
        Cl=squeeze(diffCldat.dark(ifield).fCl_arr(i,:,il));
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end  
    Cl=squeeze(diffCldat.flight(ifield).fCl_arr(:,il))';
    semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'r+','MarkerSize',3);hold off
    xlim([axlim(1),axlim(2)]);
    ylim([axlim(3),axlim(4)]);
    title(sprintf('filter'))
    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow

    subplot(2,2,3)
    for i=1:numel(dt.time)
        Cl=squeeze(diffCldat.dark(ifield).wrCl_arr(i,:,il));
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end  
    Cl=squeeze(diffCldat.flight(ifield).wrCl_arr(:,il))';
    semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'r+','MarkerSize',3);hold off
    Cl1=Cl;
    axlim=axis;
    title(sprintf('FW'))
    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    drawnow
    
    subplot(2,2,4)
    for i=1:numel(dt.time)
        Cl=squeeze(diffCldat.dark(ifield).wfCl_arr(i,:,il));
        semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'color',[0.7,0.7,0.7]);hold on
        drawnow
    end  
    Cl=squeeze(diffCldat.flight(ifield).wfCl_arr(:,il))';
    pltfilt=semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl),'r+','MarkerSize',3);hold on
    pltunfilt=semilogy(nfr_arr,sqrt(lv.*(lv+1).*Cl1),'ro','MarkerSize',5);hold off
    xlim([axlim(1),axlim(2)]);
    ylim([axlim(3),axlim(4)]);
    title(sprintf('FW + filter'))
    xlabel('$N$','interpreter','latex','fontsize',18)
    ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
        'interpreter','latex','fontsize',18)
    legend([pltfilt,pltunfilt],{'filtered','unfiltered'},'location','northeast')
    legend boxoff
    drawnow
    
    savename=sprintf('%sfield%d/ell%d',savedir,ifield,il);
    print(savename,'-dpng');close    

end
end

return 