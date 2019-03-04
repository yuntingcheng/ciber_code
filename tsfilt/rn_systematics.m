
flight=40030;
inst=1;

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/');
load(sprintf('%s/diffCldat',loaddir),'diffCldat');
%%
figure
ifield=8;
for il=9:29
l=diffCldat.l(il);

nfr_arr=diffCldat.flight(ifield).nfr_arr;
fCl_arr=diffCldat.flight(ifield).wfCl_arr(:,il)';

semilogy(nfr_arr,sqrt(l.*(l+1).*fCl_arr));hold on
drawnow
end

semilogy(nfr_arr,sqrt(l.*(l+1).*1e-9./nfr_arr./(nfr_arr.^2-1)),'k--');hold on
%% extrapolate Cl-nfr to full nfr
for ifield=4:8
figure
setwinsize(gcf,1000,600)
for il=14:21
subplot(2,4,il-13)

l=diffCldat.l(il);

nfr_arr=diffCldat.flight(ifield).nfr_arr;
nfr1_arr=2:nfr_arr(end)*2;

fCl_arr=diffCldat.flight(ifield).wfCl_arr(:,il)';
po=polyfit(nfr_arr,1./fCl_arr,2);
fClfit_arr=polyval(po,2:nfr_arr(end)*2);fClfit_arr=1./fClfit_arr;
semilogy(nfr_arr,l.*(l+1).*fCl_arr,'o','color',get_color(il));hold on
semilogy(nfr1_arr,l.*(l+1).*fClfit_arr,'color',get_color(il),'linestyle','-');

nCl_arr=mean(diffCldat.dark(ifield).wfCl_arr(:,:,il));
po=polyfit(nfr_arr,1./nCl_arr,2);
nClfit_arr=polyval(po,2:nfr_arr(end)*2);nClfit_arr=1./nClfit_arr;
semilogy(nfr_arr,l.*(l+1).*nCl_arr,'+','color',get_color(il));hold on
semilogy(nfr1_arr,l.*(l+1).*nClfit_arr,'color',get_color(il),'linestyle','--');
legend({'flgiht','flight fit','dark','dark fit'})
legend boxoff
xlabel('$N$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell$','interpreter','latex','fontsize',18)

title(sprintf('ell #%d, fCl/nCl=%.2f',il,fClfit_arr(end)./nClfit_arr(end)));
drawnow
end
savename=strcat('/Users/ytcheng/Desktop/nfrint_i',num2str(ifield));
print(savename,'-dpng');close    
end
