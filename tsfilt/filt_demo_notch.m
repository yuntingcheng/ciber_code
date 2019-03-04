%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%demo notch filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
sin_freq = 9.503;
ifield=8;
pixscale=7;

nfr=26;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',loaddir),'maskdat');
bigmask=maskdat.mask(ifield).bigmask;
%% make sin wave pick up noise
ampt_poly=[-2.5151e-5,9.5706e-4,0.0026,0.6154];
phase_in=0.71;
fr_sin=frs_sin_amp(flight,inst,nfr,sin_freq,ampt_poly,phase_in);
[raw_sin,off_sin] =linfit_map(fr_sin,'verbose',0); 
%% get signal
[~,l,~,~,lbin]=get_angular_spec(randn(1024),randn(1024),pixscale);
logx=[log10(2e3) log10(1e2)];x=10.^logx;
logy=[log10(2) log10(12)];y=10.^logy;
logCl=interp1(logx,logy,log10(l),'linear','extrap');
Clshot=ones(size(logCl)).*1e3.*2.*pi./1e5./(1e5+1);
Clin=(((10.^logCl)).*2.*pi./l./(l+1))+Clshot;
sigmap = map_from_power_spec(lbin,Clin,1024,1024,pixscale,1);
sigmap=sigmap./cal;
fr_sig=frs_from_map(sigmap,nfr,'addRN',0);
[raw_sig,off_sig] = linfit_map(fr_sig,'verbose',0);
%% get the fweight
nsim=100;
fCl2d1_stack=zeros(1024,1024);
for i=1:nsim
display(sprintf('sim%d/%d',i,nsim));

fr_read=frs_from_map(zeros(1024),nfr);
fr_tot=fr_sin+fr_read;

[raw_tot,off_tot] = linfit_map(fr_tot,'verbose',0);

[filtfr] = imager_filtts(fr_tot,off_tot,bigmask,'inst',inst,...
    'donotch',1,'dosinfilt',0,'verbose',1,'makeplot',1,...
    'notch_nus',[sin_freq],'notch_width',[0.1],'notch_iter',1); 

[filt_tot] = linfit_map(filtfr,'verbose',0);
[~,maskin1]=get_skymap(filt_tot,bigmask,4,5);

[fCl,l,~,~,binl,~,fCl2d] = get_angular_spec...
    (filt_tot.*maskin1,filt_tot.*maskin1,pixscale);

fCl2d1_stack=fCl2d1_stack+fCl2d;
end
avgfCl2d=fCl2d1_stack./nsim;
fw=(fftshift(fftshift(1./squeeze(avgfCl2d))))';
%% do the filtering
fr_read=frs_from_map(zeros(1024),nfr);
fr_tot=fr_sin+fr_read+fr_sig;

[raw_read,off_read] = linfit_map(fr_read,'verbose',0); 
[raw_tot,off_tot] = linfit_map(fr_tot,'verbose',0);

[filtfr] = imager_filtts(fr_tot,zeros(1024),bigmask,'inst',inst,...
    'donotch',1,'dosinfilt',0,'verbose',1,'makeplot',1,...
    'notch_nus',[sin_freq],'notch_width',[0.1],'notch_iter',1); 
[filt_tot] = linfit_map(filtfr,'verbose',0);

%%% filter nu=frate*512 for comparison
[filtfr1] = imager_filtts(fr_tot,zeros(1024),bigmask,'inst',inst,...
    'donotch',1,'dosinfilt',0,'verbose',1,'makeplot',1,...
    'notch_nus',[sin_freq,frate*512],'notch_width',[0.1,50],'notch_iter',1); 
[filt_tot1] = linfit_map(filtfr1,'verbose',0);
%%% sin+amp filtering for comparison
[filtfr2] = imager_filtts(fr_tot,zeros(1024),bigmask,'inst',inst,...
    'sin_freq',sin_freq,'verbose',1,'makeplot',1); 
[filt_tot2] = linfit_map(filtfr2,'verbose',0);
%% get PS
[~,maskin1]=get_skymap(filt_tot,bigmask,4,5);
[rrCl,l,~,~,binl,~,rrCl2d] = get_angular_spec...
    (raw_read.*maskin1,raw_read.*maskin1,pixscale);
[rsCl,l,~,~,binl,~,rsCl2d] = get_angular_spec...
    (raw_sig.*maskin1,raw_sig.*maskin1,pixscale);

[rCl,l,~,~,binl,~,rCl2d] = get_angular_spec...
    (raw_tot.*maskin1,raw_tot.*maskin1,pixscale);
[fCl,l,~,~,binl,~,fCl2d] = get_angular_spec...
    (filt_tot.*maskin1,filt_tot.*maskin1,pixscale);
[fwCl,l,~,~,binl,~,fwCl2d] = get_angular_spec...
    (filt_tot.*maskin1,filt_tot.*maskin1,pixscale,'w',fw);

[~,maskin1]=get_skymap(filt_tot1,bigmask,4,5);
[fCl1,l,~,~,binl,~,fCl2d1] = get_angular_spec...
    (filt_tot1.*maskin1,filt_tot1.*maskin1,pixscale);
[~,maskin1]=get_skymap(filt_tot2,bigmask,4,5);
[fCl2,l,~,~,binl,~,fCl2d2] = get_angular_spec...
    (filt_tot2.*maskin1,filt_tot2.*maskin1,pixscale);
%% comparing sin+amp vs notch 2D
ell=get_l(1024,1024,pixscale,1);
fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
[x,y]=find(fmask);

figure
setwinsize(gcf,1200,1200)

subplot(3,3,1)
imageclip(raw_tot);
v=caxis;
title('unfilt map');

subplot(3,3,2)
imageclip(filt_tot2);
caxis(v);
title('sin+amp filter');

subplot(3,3,3)
imageclip(filt_tot);
caxis(v);
title('notch filter');

subplot(3,3,4)
imageclip(rCl2d);
v=caxis;
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('unfilt 2D Cl');

subplot(3,3,5)
imageclip(fCl2d2);
caxis(v);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('sin+amp filter 2D Cl');

subplot(3,3,6)
imageclip(fCl2d);
caxis(v);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('notch filter 2D Cl');


errCl2d=(rCl2d-(rsCl2d+rrCl2d))./(rsCl2d+rrCl2d);
subplot(3,3,7)
imageclip(errCl2d);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('[unfilt-(signal+white RN)]/[signal+white RN]')
caxis([-2,2])

errCl2d=(fCl2d2-(rsCl2d+rrCl2d))./(rsCl2d+rrCl2d);
subplot(3,3,8)
imageclip(errCl2d);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('[sinfilt-(signal+white RN)]/[signal+white RN]')
caxis([-2,2])

errCl2d=(fCl2d-(rsCl2d+rrCl2d))./(rsCl2d+rrCl2d);
subplot(3,3,9)
imageclip(errCl2d);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('[notchfilt-(signal+white RN)]/[signal+white RN]')
caxis([-2,2])
%% comparing sin+amp vs notch 1D
figure
loglog(l,l.*(l+1).*rCl,'-','DisplayName','unfiltered:signal+ white RN');hold on
loglog(l,l.*(l+1).*fCl2,'bo','DisplayName','sin filter');hold on
loglog(l,l.*(l+1).*fCl,'ro','DisplayName','notch filter');hold on
loglog(l,l.*(l+1).*(rsCl+rrCl),'k','DisplayName','signal+white RN');hold on
ylim([1e-4,1e-1]);
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell$','interpreter','latex','fontsize',18)
%% plot 2DPS notch
ell=get_l(1024,1024,pixscale,1);
fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
[x,y]=find(fmask);

figure
setwinsize(gcf,1500,800)

subplot(2,3,1)
imageclip(raw_tot);
v=caxis;
title('unfilt map');

subplot(2,3,2)
imageclip(filt_tot);
caxis(v);
title('filt map');

subplot(2,3,4)
imageclip(rCl2d);
v=caxis;
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('unfilt 2D Cl');

subplot(2,3,5)
imageclip(fCl2d);
caxis(v);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('filt 2D Cl');

drawnow

errCl2d=(fCl2d-(rsCl2d+rrCl2d))./(rsCl2d+rrCl2d);
subplot(2,3,3)
imageclip(errCl2d);
title('[filt-(signal+white RN)]/[signal+white RN]')
caxis([-2,2])

subplot(2,3,6);
imageclip(errCl2d(min(x):max(x),min(y):max(y)));hold on
caxis([-2,2])
[c,h]=contour(ell(min(x):max(x),min(y):max(y)),[binl(20),binl(21)],'color','m');
set(h,'linewidth',1.5);
[c,h]=contour(ell(min(x):max(x),min(y):max(y)),...
    [binl(18),binl(19),binl(22),binl(23),binl(24)],'color','b');
set(h,'linewidth',1.5);
title('[filt-(signal+white RN)]/[signal+white RN]')
%% plot 1DPS notch
figure
loglog(l,l.*(l+1).*rCl,'-','DisplayName','unfiltered:signal+ white RN');hold on
loglog(l,l.*(l+1).*fCl,'ro','DisplayName','filter 9.503 Hz');hold on
loglog(l,l.*(l+1).*fCl1,'b+','DisplayName','filter 9.503+287 Hz');hold on
loglog(l,l.*(l+1).*fwCl,'r+','DisplayName','filtered+ FW');hold on
loglog(l,l.*(l+1).*(rsCl+rrCl),'k-','DisplayName','signal+white RN');hold on
loglog(l,l.*(l+1).*(rsCl),'linestyle','--','color',[0,.8,.2],'DisplayName','signal');
loglog(l,l.*(l+1).*(rrCl),'linestyle','-.','color',[0,.8,.2],'DisplayName','white RN');
ylim([1e-4,1e-1]);
xlim([1e2,2e5]);
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell$','interpreter','latex','fontsize',18)
