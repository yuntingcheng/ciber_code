%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot 25 frames auto PS from get_auto25ps.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band=2;

psdatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(psdatdir,'band',num2str(band),'_fullpsdat'),'psdat');
rndir='/Users/ytcheng/ciber/doc/20160905_NoiseModel/full/PS1D/';
phdir=strcat('/Users/ytcheng/ciber/doc/20160919_FlightDat/band',...
    num2str(band),'_ph/');
savedir=strcat('/Users/ytcheng/ciber/doc/20160919_FlightDat/band',...
    num2str(band),'_fullps/');
%%
for i=1:8
l=psdat(i).l;
fig=figure;

%%% raw
yvalue=(l.*(l+1).*psdat(i).rCl./2./pi);
ebar=(l.*(l+1).*psdat(i).drCl./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltrCl=errorbar(l,yvalue,ebarlow,ebar,'.','color',[0.6,0.6,0.6]...
                            ,'markersize',20);hold on


%%% before RN sub
yvalue=(l.*(l+1).*psdat(i).sCl5./2./pi);
ebar=(l.*(l+1).*psdat(i).dsCl5./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltCl=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',20);hold on

%%% after RN sub
yvalue=(l.*(l+1).*psdat(i).subsCl5./2./pi);
ebar=(l.*(l+1).*psdat(i).subdsCl5./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltsubCl=errorbar(l,yvalue,ebarlow,ebar,'.r','markersize',20);hold on

%%% read noise
load(sprintf('%sb%d_i%d_wnClFclip_arr',rndir,band,i),'wnClFclip_arr');

y1=(l.*(l+1).*(prctile(wnClFclip_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(wnClFclip_arr,84))./2./pi);
y1=y1(9:29);y2=y2(9:29);ll=l(9:29);
pltrn=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0,0,1],'facealpha',0.2,'EdgeColor','none');hold on

y1=(l.*(l+1).*min(wnClFclip_arr)./2./pi);
y2=(l.*(l+1).*max(wnClFclip_arr)./2./pi);
y1=y1(9:29);y2=y2(9:29);ll=l(9:29);
pltrnmm=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0,0.8,1],'facealpha',0.2,'EdgeColor','none');hold on

%%% photon noise
load(strcat(phdir,'i',num2str(i),'phCl_arr'),'phCl_arr');
y1=(l.*(l+1).*(prctile(phCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(phCl_arr,84))./2./pi);
y1=y1(9:29);y2=y2(9:29);ll=l(9:29);
pltph=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',1,'EdgeColor','none');


xlim([1e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
title(strcat(psdat(i).name));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltrCl,pltCl,pltsubCl,pltrn,pltph],{'raw','cal',...
    'cal-RN','RN','photon noise'},...
       'Location','southeast','FontSize',15);
legend boxoff
imname=strcat(savedir,'b',num2str(band),'_i',num2str(i),'_fullps');
print(imname,'-dpng');%close
end
%% simplified version
for i=1:8
l=psdat(i).l;
fig=figure;

%%% read noise + photon noise
load(sprintf('%sb%d_i%d_wnClFclip_arr',rndir,band,i),'wnClFclip_arr');
load(strcat(phdir,'i',num2str(i),'phCl_arr'),'phCl_arr');

y1=(l.*(l+1).*(prctile(wnClFclip_arr,8.8));
y2=(l.*(l+1).*(prctile(wnClFclip_arr,91.2));
y1=y1(9:29);y2=y2(9:29);ll=l(9:29);
pltn=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

%%% before RN sub
yvalue=(l.*(l+1).*psdat(i).sCl5./2./pi);
ebar=(l.*(l+1).*psdat(i).dsCl5./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltCl=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',20);hold on

xlim([1e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
title(strcat(psdat(i).name));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltCl,pltn],{'cal','noise'},...
       'Location','southeast','FontSize',15);
legend boxoff
imname=strcat(savedir,'b',num2str(band),'_i',num2str(i),'_fullps_simp');
print(imname,'-dpng');%close
end
%% no Fourier weighting (simplified version)
for i=1:8
l=psdat(i).l;
fig=figure;

%%% read noise + photon noise
load(sprintf('%sb%d_i%d_nClF_arr',rndir,band,i),'nClF_arr');
load(strcat(phdir,'i',num2str(i),'phCl_arr'),'phCl_arr');

y1=(l.*(l+1).*(prctile(nClF_arr,8.8)+prctile(phCl_arr,8.8))./2./pi);
y2=(l.*(l+1).*(prctile(nClF_arr,91.2)+prctile(phCl_arr,91.2))./2./pi);
y1=y1(9:29);y2=y2(9:29);ll=l(9:29);
pltn=fill([ll,flip(ll)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

%%% before RN sub
yvalue=(l.*(l+1).*psdat(i).sCl5nw./2./pi);
ebar=(l.*(l+1).*psdat(i).dsCl5nw./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltCl=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',20);hold on

xlim([1e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
title(strcat(psdat(i).name));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltCl,pltn],{'cal','noise'},...
       'Location','southeast','FontSize',15);
legend boxoff
imname=strcat(savedir,'b',num2str(band),'_i',num2str(i),'_fullps_simp_nw');
print(imname,'-dpng');%close
end
