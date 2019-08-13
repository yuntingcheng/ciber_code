flight=40030;
inst=2;
ifield=5;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

m_min = 19;
m_max = 20;
%% Comapring stacking all sources / uniform sources
figure
setwinsize(gcf,800,600)

sample_type = 'jack_random';
load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s',...
    savedir,dt.name,m_min,m_max,sample_type),'stackdat');

subplot(2,2,1)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profcbs,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profcbg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profcbs,stackdat.errjack.profcbs,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profcbg,stackdat.errjack.profcbg,'b');
xlim([4e-1,1.1e3])
ylim([-1,10])
title('CIBER, all srcs','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

subplot(2,2,3)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profpss,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profpsg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profpss,stackdat.errjack.profpss,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profpsg,stackdat.errjack.profpsg,'b');
xlim([4e-1,1.1e3])
ylim([-0.05,0.1])
title('PanSTARRS, all srcs','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

sample_type = 'jack_random_uni';
load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s',...
    savedir,dt.name,m_min,m_max,sample_type),'stackdat');

subplot(2,2,2)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profcbs,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profcbg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profcbs,stackdat.errjack.profcbs,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profcbg,stackdat.errjack.profcbg,'b');
xlim([4e-1,1.1e3])
ylim([-1,10])
title('CIBER, uniform srcs','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

subplot(2,2,4)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profpss,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profpsg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profpss,stackdat.errjack.profpss,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profpsg,stackdat.errjack.profpsg,'b');
xlim([4e-1,1.1e3])
ylim([-0.05,0.1])
title('PanSTARRS, uniform srcs','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

s=suptitle(sprintf('TM%d %s %d<m<%d',inst,dt.name,m_min,m_max))
set(s,'FontSize',15);
%% Comparing error from random and regional jackknife
figure
setwinsize(gcf,1200,600)

sample_type = 'jack_random';
load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s',...
    savedir,dt.name,m_min,m_max,sample_type),'stackdat');
cbgerr_rand = stackdat.errjack.profcbg;
psgerr_rand = stackdat.errjack.profpsg;
cbserr_rand = stackdat.errjack.profcbs;
psserr_rand = stackdat.errjack.profpss;

subplot(2,3,1)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profcbs,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profcbg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profcbs,stackdat.errjack.profcbs,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profcbg,stackdat.errjack.profcbg,'b');
xlim([4e-1,1.1e3])
ylim([-1,10])
title('CIBER, all srcs, random jackknife','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

subplot(2,3,4)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profpss,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profpsg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profpss,stackdat.errjack.profpss,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profpsg,stackdat.errjack.profpsg,'b');
xlim([4e-1,1.1e3])
ylim([-0.05,0.1])
title('PanSTARRS, all srcs, random jackknife','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

sample_type = 'jack_region';
load(sprintf('%s/stackdat_%s_%d_%d_%s',...
    savedir,dt.name,m_min,m_max,sample_type),'stackdat');
cbgerr_region = stackdat.errjack.profcbg;
psgerr_region = stackdat.errjack.profpsg;
cbserr_region = stackdat.errjack.profcbs;
psserr_region = stackdat.errjack.profpss;

subplot(2,3,2)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profcbs,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profcbg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profcbs,stackdat.errjack.profcbs,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profcbg,stackdat.errjack.profcbg,'b');
xlim([4e-1,1.1e3])
ylim([-1,10])
title('CIBER, all srcs, subregion jackknife','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

subplot(2,3,5)
semilogx(stackdat.r_arr.*0.99,stackdat.all.profpss,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profpsg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profpss,stackdat.errjack.profpss,'r');
errorbar(stackdat.r_arr,stackdat.all.profpsg,stackdat.errjack.profpsg,'b');
xlim([4e-1,1.1e3])
ylim([-0.05,0.1])
title('PanSTARRS, all srcs, subregion jackknife','fontsize',15);
xlabel('arcsec', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)


subplot(2,3,3)
semilogx(stackdat.r_arr.*0.99,cbserr_region./cbserr_rand,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,cbgerr_region./cbgerr_rand,'b.-');hold on
legend({'stars','galaxies'});
xlim([4e-1,1.1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('subregion err / random err', 'fontsize',15);
title('error ratio', 'fontsize',15);

subplot(2,3,6)
semilogx(stackdat.r_arr.*0.99,psserr_region./psserr_rand,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,psgerr_region./psgerr_rand,'b.-');hold on
legend({'stars','galaxies'});
xlim([4e-1,1.1e3])
xlabel('arcsec', 'fontsize',15);
ylabel('subregion err / random err', 'fontsize',15);
title('error ratio','fontsize',15);

s=suptitle(sprintf('TM%d %s %d<m<%d',inst,dt.name,m_min,m_max));
set(s,'FontSize',15);
%% plot CIBER / PanSTARRS smoothed image

flight=40030;
inst=1;
ifiel=4;

mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
dx = 1200;
verbose = 0;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;
m_min_arr = m_min_arr(10:13);
m_max_arr = m_max_arr(10:13);

figure
setwinsize(gcf,1100,600)
smscale = [10,50,100,200];

subplot('Position',[0.03,0.65,0.16,0.3])
imageclip(cbmap.*mask_inst.*strmask);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
ylabel('CIBER','fontsize',15);
title('no smoothing','fontsize',15);
subplot('Position',[0.03,0.35,0.16,0.3])
imageclip(psmap.*mask_inst.*strmask);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
ylabel('PanSTAARS','fontsize',15);
subplot('Position',[0.03,0.05,0.16,0.3])
imageclip((cbmap-psmap).*mask_inst.*strmask);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
ylabel('CIBER - PanSTARRS','fontsize',15);

for i=1:numel(smscale)
    sm = smscale(i);
    smcb = fillpadsmooth(cbmap.*mask_inst.*strmask,mask_inst.*strmask,sm);
    smps = fillpadsmooth(psmap.*mask_inst.*strmask,mask_inst.*strmask,sm);
  
    subplot('Position',[0.03+i*0.2,0.65,0.16,0.3])
    imageclip(smcb);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    title(sprintf('%d{''}{''} Gaussian smoothing',sm*7),'fontsize',15);
    subplot('Position',[0.03+i*0.2,0.35,0.16,0.3])
    imageclip(smps);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    subplot('Position',[0.03+i*0.2,0.05,0.16,0.3])
    imageclip(smcb-smps);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
end
%% plot HSC smoothed image
flight=40030;
inst=1;
ifield=8;
hsc_idx=2;
f_ihl=0;
mypaths=get_paths(flight);
name = HSC_fields_info(hsc_idx);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackmapdathsc_%s%d',loaddir,name,f_ihl*100),...
            'stackmapdatsim');
cbmap = stackmapdatsim(ifield).all.cbmap;
psmap = stackmapdatsim(ifield).all.psmap;
mask_inst = stackmapdatsim(ifield).all.mask_inst_clip;
strmask = stackmapdatsim(ifield).all.strmask;

figure
setwinsize(gcf,1100,600)
smscale = [10,50,100,200];

subplot('Position',[0.03,0.65,0.16,0.3])
imageclip(cbmap.*mask_inst.*strmask);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
ylabel('HSC all','fontsize',15);
title('no smoothing','fontsize',15);
subplot('Position',[0.03,0.35,0.16,0.3])
imageclip(psmap.*mask_inst.*strmask);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
ylabel('HSC m<20','fontsize',15);
subplot('Position',[0.03,0.05,0.16,0.3])
imageclip((cbmap-psmap).*mask_inst.*strmask);
set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
ylabel('HSC all - HSC m<20','fontsize',15);

for i=1:numel(smscale)
    sm = smscale(i);
    smcb = fillpadsmooth(cbmap.*mask_inst.*strmask,mask_inst.*strmask,sm);
    smps = fillpadsmooth(psmap.*mask_inst.*strmask,mask_inst.*strmask,sm);
  
    subplot('Position',[0.03+i*0.2,0.65,0.16,0.3])
    imageclip(smcb);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    title(sprintf('%d{''}{''} Gaussian smoothing',sm*7),'fontsize',15);
    subplot('Position',[0.03+i*0.2,0.35,0.16,0.3])
    imageclip(smps);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
    subplot('Position',[0.03+i*0.2,0.05,0.16,0.3])
    imageclip(smcb-smps);
    set(gca,'xtick',[],'xticklabel',[],'ytick',[],'yticklabel',[])
end
%% plot profile with smoothing
flight=40030;
inst=1;
ifield=4;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

m_min = 19;
m_max = 20;

figure
setwinsize(gcf,1100,400)

sample_type = 'jack_random';
load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s',...
    savedir,dt.name,m_min,m_max,sample_type),'stackdat');

subplot('Position',[0.04,0.60,0.16,0.35])
semilogx(stackdat.r_arr.*0.99,stackdat.all.profcbs,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profcbg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profcbs,stackdat.errjack.profcbs,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profcbg,stackdat.errjack.profcbg,'b');
xlim([4e-1,1.1e3])
ylim([-1,10])
title('CIBER','fontsize',10);
xlabel('arcsec', 'fontsize',10)
ylabel('I [nW/m^2/sr]', 'fontsize',10)

subplot('Position',[0.04,0.10,0.16,0.35])
semilogx(stackdat.r_arr.*0.99,stackdat.all.profpss,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profpsg,'b.-');hold on
legend({'stars','galaxies'});
errorbar(stackdat.r_arr.*0.99,stackdat.all.profpss,stackdat.errjack.profpss,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profpsg,stackdat.errjack.profpsg,'b');
xlim([4e-1,1.1e3])
ylim([-0.05,0.1])
title('PanSTARRS','fontsize',10);
xlabel('arcsec', 'fontsize',10)
ylabel('I [nW/m^2/sr]', 'fontsize',10)

ism=0;
for sm_scale = [10,50,100,200]
ism=ism+1;
load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s_sm%d',...
    savedir,dt.name,m_min,m_max,sample_type,sm_scale),'stackdat');

subplot('Position',[0.04+ism*0.19,0.60,0.16,0.35])
semilogx(stackdat.r_arr.*0.99,stackdat.all.profcbs,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profcbg,'b.-');hold on
legend({'stars','galaxies'});
vline(sm_scale*7,'k--');
errorbar(stackdat.r_arr.*0.99,stackdat.all.profcbs,stackdat.errjack.profcbs,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profcbg,stackdat.errjack.profcbg,'b');
xlim([4e-1,1.1e3])
ylim([-1,10])
title(sprintf('CIBER, smooth %d pix',sm_scale),'fontsize',10);
xlabel('arcsec', 'fontsize',10)

subplot('Position',[0.04+ism*0.19,0.10,0.16,0.35])
semilogx(stackdat.r_arr.*0.99,stackdat.all.profpss,'r.-');hold on
semilogx(stackdat.r_arr.*1.01,stackdat.all.profpsg,'b.-');hold on
legend({'stars','galaxies'});
vline(sm_scale*7,'k--');
errorbar(stackdat.r_arr.*0.99,stackdat.all.profpss,stackdat.errjack.profpss,'r');
errorbar(stackdat.r_arr.*1.01,stackdat.all.profpsg,stackdat.errjack.profpsg,'b');
xlim([4e-1,1.1e3])
ylim([-0.05,0.1])
title(sprintf('PanSTARRS, smooth %d pix',sm_scale),'fontsize',10);
xlabel('arcsec', 'fontsize',10)

end
s=suptitle(sprintf('TM%d %s %d<m<%d',inst,dt.name,m_min,m_max));
set(s,'FontSize',15);
%% plot exceed
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
sm_scale_arr = [-1,0,10,50,100,200];
for ism=1:6
figure
setwinsize(gcf,1400,300)

    sm_scale = sm_scale_arr(ism);
for im=1:4
    m_min = im + 15;
    m_max = m_min + 1;
    subplot(1,4,im)
    name = {};
    vavg = 0;
    eavg = 0;
    for ifield = 4:8
        dt=get_dark_times(flight,inst,ifield);
        name{end+1} = dt.name;
        if sm_scale==0
            load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s',...
                savedir,dt.name,m_min,m_max,sample_type),'stackdat');
        elseif sm_scale==-1
            load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s',...
                savedir,dt.name,m_min,m_max,'jack_random_uni'),'stackdat');
        else
            load(sprintf('%s/TM1_tests/stackdat_%s_%d_%d_%s_sm%d',...
                savedir,dt.name,m_min,m_max,sample_type,sm_scale),'stackdat');
        end
        
        v = (stackdat.all.profcbg100 - stackdat.all.profcbs100)...
            - (stackdat.all.profpsg100 - stackdat.all.profpss100);
        e = sqrt(stackdat.errjack.profcbs100^2+stackdat.errjack.profcbg100^2 ...
            +stackdat.errjack.profpss100^2+stackdat.errjack.profpsg100^2);
        errorbar(ifield-3,v,e,'k.','markersize',10);hold on
        vavg = vavg + v/(e^2);
        eavg = eavg + 1/(e^2);
        
    end
    vavg = vavg/eavg;
    eavg = sqrt(1/eavg);
    legend({sprintf('average: %.3f +- %.3f \n SNR = %.3f',vavg,eavg,abs(vavg/eavg))})
    plot([0,6],[0,0],'b--')
    ylabel('mean excess > 100 arcsec [nW/m^2/sr]', 'fontsize',10);
    xlim([0.5,5.5])
    ylim([-0.5,0.8])
    set(gca,'xTick',1:5);
    set(gca, 'xTickLabels', name); 
    title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);
    

end
if sm_scale>0
    s=suptitle(sprintf('smooth %d pix',sm_scale));
    set(s,'FontSize',15);
elseif sm_scale==-1
    s=suptitle('uniform sources');
    set(s,'FontSize',15);
end

print(sprintf('/Users/ytcheng/Desktop/excess_sm%d',sm_scale),'-dpng');
end