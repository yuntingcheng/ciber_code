function plot_Cl_l(flight,inst)
mypaths=get_paths(flight);

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/diffplot/');
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/');
load(sprintf('%s/diffCldat',loaddir),'diffCldat');

l=diffCldat.l;
for ifield=1:8

dt=get_dark_times(flight,inst,ifield);
nfr_arr=diffCldat.dark(ifield).nfr_arr;

for infr=1:numel(nfr_arr)
figure
setwinsize(gcf,1000,1000)
subplot(2,2,1)
for i=1:numel(dt.time)
    Cl=squeeze(diffCldat.dark(ifield).rCl_arr(i,infr,:))';
    if i~=1
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    else
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    end        
end
Cl=squeeze(diffCldat.flight(ifield).rCl_arr(infr,:));
loglog(l,sqrt(l.*(l+1).*Cl),'r+','MarkerSize',3);hold off 
xlim([1e2,2e5]);
axlim=axis;
title(sprintf('%s, nfr=%d',dt.name,nfr_arr(infr)))
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
    'interpreter','latex','fontsize',18)
drawnow

subplot(2,2,2)
for i=1:numel(dt.time)
    Cl=squeeze(diffCldat.dark(ifield).fCl_arr(i,infr,:))';
    if i~=1
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    else
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    end        
end
Cl=squeeze(diffCldat.flight(ifield).fCl_arr(infr,:));
loglog(l,sqrt(l.*(l+1).*Cl),'r+','MarkerSize',3);hold off 
xlim([axlim(1),axlim(2)]);
ylim([axlim(3),axlim(4)]);
title(sprintf('filter'))
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
    'interpreter','latex','fontsize',18)
drawnow

subplot(2,2,3)
for i=1:numel(dt.time)
    Cl=squeeze(diffCldat.dark(ifield).wrCl_arr(i,infr,:))';
    if i~=1
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    else
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    end        
end
Cl=squeeze(diffCldat.flight(ifield).wrCl_arr(infr,:));
Cl1=Cl;
loglog(l,sqrt(l.*(l+1).*Cl),'r+','MarkerSize',3);hold off 
xlim([axlim(1),axlim(2)]);
ylim([axlim(3),axlim(4)]);
title(sprintf('FW'))
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
    'interpreter','latex','fontsize',18)
drawnow

subplot(2,2,4)
for i=1:numel(dt.time)
    Cl=squeeze(diffCldat.dark(ifield).wfCl_arr(i,infr,:))';
    if i~=1
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    else
    loglog(l,sqrt(l.*(l+1).*Cl),'color',[0.7,0.7,0.7]);hold on
    drawnow
    end        
end
Cl=squeeze(diffCldat.flight(ifield).wfCl_arr(infr,:));
pltfilt=loglog(l,sqrt(l.*(l+1).*Cl),'r+','MarkerSize',3);hold on
pltunfilt=loglog(l,sqrt(l.*(l+1).*Cl1),'ro','MarkerSize',5);hold off 

xlim([axlim(1),axlim(2)]);
ylim([axlim(3),axlim(4)]);
title(sprintf('FW + filter'))
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\sqrt{\ell(\ell+1)C_\ell}$',...
    'interpreter','latex','fontsize',18)
legend([pltfilt,pltunfilt],{'filtered','unfiltered'},'location','southeast')
legend boxoff
drawnow
savename=sprintf('%sfield%d/nfr%d',savedir,ifield,nfr_arr(infr));
print(savename,'-dpng');close    

end
end

return 