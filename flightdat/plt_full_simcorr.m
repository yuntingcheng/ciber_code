%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot 25 frames auto PS from get_auto25ps.m, 
%corrected with simulation scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
beam=6.2;
band=1;

psdatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(psdatdir,'band',num2str(band),'_fullpsdat'),'psdat');
rndir='/Users/ytcheng/ciber/doc/20160905_NoiseModel/full/PS1D/';
phdir=strcat('/Users/ytcheng/ciber/doc/20160919_FlightDat/band',...
    num2str(band),'_ph/');
savedir=strcat('/Users/ytcheng/ciber/doc/20160919_FlightDat/band',...
    num2str(band),'_fullps_simc/');
simdir='/Users/ytcheng/ciber/doc/20160921_Simulation/';
%%
for i=1:8
load(strcat(simdir,'Cldat/i',num2str(i),'_rClin_arr'),'rClin_arr');

l=psdat(i).l;
fig=figure;

%%% before RN sub
yvalue=(l.*(l+1).*psdat(i).sCl5./2./pi);
ebarup=(l.*(l+1).*psdat(i).dsCl5./2./pi);
ebarlow=ebarup;
pltCl=errorbar(l,yvalue,ebarlow,ebarup,'.k','markersize',20);hold on

%%% rescaling
yvalue=(l.*(l+1).*(psdat(i).sCl5./prctile(rClin_arr,50))./2./pi);
ebarup=yvalue-(l.*(l+1).*(psdat(i).sCl5./prctile(rClin_arr,84))./2./pi);
ebarlow=(l.*(l+1).*(psdat(i).sCl5./prctile(rClin_arr,16))./2./pi)-yvalue;
ebarlow(find(ebarlow>yvalue))=yvalue(find(ebarlow>yvalue)).*(1-1e-10);
pltClc=errorbar(l,yvalue,ebarlow,ebarup,'.r','markersize',20);hold on

%%% read noise
load(strcat(rndir,'b',num2str(band),...
        '_i',num2str(i),'_nClfull_arr'),'nClfull_arr');

y1=(l.*(l+1).*(prctile(nClfull_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(nClfull_arr,84))./2./pi);
y1=y1(10:29);y2=y2(10:29);ll=l(10:29);
pltrn=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0,0,1],'facealpha',0.2,'EdgeColor','none');hold on

y1=(l.*(l+1).*min(nClfull_arr)./2./pi);
y2=(l.*(l+1).*max(nClfull_arr)./2./pi);
y1=y1(10:29);y2=y2(10:29);ll=l(10:29);
pltrnmm=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0,0.8,1],'facealpha',0.2,'EdgeColor','none');hold on

%%% photon noise
load(strcat(phdir,'i',num2str(i),'phCl_arr'),'phCl_arr');
y1=(l.*(l+1).*(prctile(phCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(phCl_arr,84))./2./pi);
y1=y1(10:29);y2=y2(10:29);ll=l(10:29);
pltph=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',1,'EdgeColor','none');


xlim([2e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
title(strcat(psdat(i).name));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltCl,pltClc,pltrn,pltph],{'cal','sim-corr','RN','phN'},...
       'Location','southeast','FontSize',15);
legend boxoff
imname=strcat(savedir,'b',num2str(band),'_i',num2str(i),'_fullps_simc');
print(imname,'-dpng');close
end
