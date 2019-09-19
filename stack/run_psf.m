%%%% run bright source stack for a range of mag bins %%%%%%

flight = 40030;
inst = 2;

mypaths=get_paths(flight);
load(sprintf('%s/TM%d/stackmapdat',mypaths.alldat,1),'stackmapdat');
stackmapdat1 = stackmapdat;
load(sprintf('%s/TM%d/stackmapdat',mypaths.alldat,2),'stackmapdat');
stackmapdat2 = stackmapdat;

if inst==1
    stackmapdat = stackmapdat1;
else
    stackmapdat = stackmapdat2;
end
  
dx = 1200;
verbose = false;
m_min_arr = [4,4,4,4,12,13,15];
m_max_arr = [9,10,11,12,13,14,16];

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
   
for im= 1:numel(m_min_arr)

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    mask_inst = zeros([2,1024,1024]);
    mask_inst(1,:,:) = stackmapdat1(ifield).mask_inst_clip;
    mask_inst(2,:,:) = stackmapdat2(ifield).mask_inst_clip;
    strmask = stackmapdat(ifield).strmask;
    strnum = stackmapdat(ifield).strnum;    
    
    psfdat.m_min = m_min;
    psfdat.m_max = m_max;
    
    srcdat = tm_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
    'sample_type','jack_random','Nsub',10);
    
    if srcdat.Ns==0
        continue
    end

    [clipmaxs, clipmins, r_arr]=...
    stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
    mask_inst,strnum,1000,verbose,[],nan,true);
    psfdat.r_arr = r_arr;
    mask_inst = squeeze(mask_inst(inst,:,:));
    
    cbmean = mean(cbmap(find(mask_inst.*strmask)));
    psmean = mean(psmap(find(mask_inst.*strmask)));

    for isub=1:10
        [~,~,~,profcbs,profpss,profhits] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,nan,clipmaxs,clipmins,...
            srcdat.sub(isub).xs_arr,srcdat.sub(isub).ys_arr,...
            srcdat.sub(isub).ms_arr,true);

        fprintf('stack %s, %d<m<%d, isub %d, %d srcs\n',...
            dt.name,m_min,m_max,isub,srcdat.sub(isub).Ns);

        psfdat.sub(isub).counts = srcdat.sub(isub).Ns;
        profcbs(profhits==0) = 0;
        profpss(profhits==0) = 0;
        psfdat.sub(isub).profcbs = profcbs - cbmean;
        psfdat.sub(isub).profpss = profpss - psmean;        
        psfdat.sub(isub).profhits = profhits;
    end
    
    %%% profile combining all subset
    profcbs = zeros(size(psfdat.r_arr));
    profpss = zeros(size(psfdat.r_arr));
    profhits = zeros(size(psfdat.r_arr));
    counts = 0;
    for isub=1:10
        profcbs = profcbs + ...
            psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        profpss = profpss + ...
            psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        profhits = profhits + psfdat.sub(isub).profhits;
        counts = counts + psfdat.sub(isub).counts;
    end
    psfdat.all.profcbs = profcbs./profhits;
    psfdat.all.profpss = profpss./profhits;
    psfdat.all.counts = counts;
    
    %%% profile of jackknife samples (leave one out)
    for isub=1:10
        jackcbs = profcbs - psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        jackpss = profpss - psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        jackhits = profhits - psfdat.sub(isub).profhits;
        psfdat.jack(isub).profcbs = jackcbs./jackhits;
        psfdat.jack(isub).profpss = jackpss./jackhits;
    end
    
    %%% error bar with jackknife
    errcbs = zeros(size(psfdat.r_arr));
    errpss = zeros(size(psfdat.r_arr));
    for isub=1:10
        errcbs = errcbs + ...
            (psfdat.jack(isub).profcbs - psfdat.all.profcbs).^2;        
        errpss = errpss + ...
            (psfdat.jack(isub).profpss - psfdat.all.profpss).^2;
    end
    psfdat.errjack.profcbs = sqrt(errcbs.*(9/10));
    psfdat.errjack.profpss = sqrt(errpss.*(9/10));
    
    psfdatall(im).psfdat = psfdat;
end
psfdatallfields(ifield).psfdatall = psfdatall;
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/psfdat',savedir),'psfdatallfields');

%% plot the stack and get the combined PSF
flight = 40030;
mypaths=get_paths(flight);
m_min_arr = [4,4,4,4,12,13,15];
m_max_arr = [9,10,11,12,13,14,16];

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);

for inst=[1,2]
figure
setwinsize(gcf, 1500, 300)
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/psfdat',savedir),'psfdatallfields');
psfdatall = psfdatallfields(ifield).psfdatall;
for im= 1:numel(m_min_arr)
    off = 0.92 + im*0.02;
    psfdat = psfdatall(im).psfdat;
    r_arr = psfdat.r_arr;
    subplot(1,3,1)
    loglog(r_arr.*off, psfdat.all.profcbs,'.-','color',get_color(im),...
    'DisplayName',sprintf('%d < m_J < %d (%d sources)', ...
    psfdat.m_min, psfdat.m_max, psfdat.all.counts));hold on
end
vline(r_arr(11),'k--');
h=legend('show','Location','southwest');
set(h,'fontsize',8)
legend boxoff
grid on
xlim([4e-1,1.1e3])
ylim([1e-2,1.1e6])
xlabel('r [arcsec]', 'fontsize',15)
ylabel('I [nW/m^2/sr]', 'fontsize',15)

snrs = [];
for im= 1:numel(m_min_arr)
    off = 0.92 + im*0.02;
    subplot(1,3,1)
    psfdat = psfdatall(im).psfdat;
    loglog(r_arr.*off, -psfdat.all.profcbs,'o','color',get_color(im));
    errorbar(r_arr.*off, psfdat.all.profcbs, psfdat.errjack.profcbs, ...
        '.','color',get_color(im));
    errorbar(r_arr.*off, -psfdat.all.profcbs, psfdat.errjack.profcbs,...
        'o','color',get_color(im));
    if im<=4
        subplot(1,3,2)
        snr = psfdat.all.profcbs./psfdat.errjack.profcbs;
        snrs = [snrs; snr(15:18)];
        semilogx(r_arr.*off, snr, ...
            '.-','color',get_color(im),'markersize',10);hold on
        ylim([-1,15])
        xlim([4e-1,1.1e3])
        grid on
        xlabel('r [arcsec]', 'fontsize',15)
        ylabel('SNR', 'fontsize',15)
    end    
end
subplot(1,3,2)
bandname = ['I','H'];
title(sprintf('%s band,%s',bandname(inst),dt.name),'fontsize',15);

subplot(1,3,3)

im = 6;
psfdat = psfdatall(im).psfdat;
norm = psfdat.all.profcbs(1);
psfin = psfdat.all.profcbs(1:11)./norm;
psfin_err = psfdat.errjack.profcbs(1:11)./norm;
normps = psfdat.all.profpss(1);
psfinps = psfdat.all.profpss(1:11)./normps;
psfinps_err = psfdat.errjack.profpss(1:11)./normps;
loglog(r_arr(1:11), psfin,'.-','color',get_color(im));hold on
loglog(r_arr(1:11), -psfin,'o','color',get_color(im));
errorbar(r_arr(1:11), psfin, psfin_err, '.','color',get_color(im));
errorbar(r_arr(1:11), -psfin, psfin_err,'o','color',get_color(im));

psfdatallfields(ifield).r_arr = r_arr;
psfdatallfields(ifield).idx_r_cut = 11;
psfdatallfields(ifield).r_cut = r_arr(11);
psfdatallfields(ifield).im_out = im;
psfdatallfields(ifield).m_min_out = psfdat.m_min;
psfdatallfields(ifield).m_max_out = psfdat.m_max;

[~,im] = max(sum(snrs'));
norm  = norm./psfdat.all.profcbs(11);
normps  = normps./psfdat.all.profpss(11);
psfdat = psfdatall(im).psfdat;
norm = norm.*psfdat.all.profcbs(11);
psfout = psfdat.all.profcbs(12:end)./norm;
psfout_err = psfdat.errjack.profcbs(12:end)./norm;
normps = normps.*psfdat.all.profpss(11);
psfoutps = psfdat.all.profpss(12:end)./normps;
psfoutps_err = psfdat.errjack.profpss(12:end)./normps;

loglog(r_arr(11:end), [psfin(end),psfout],'.-','color',get_color(im));hold on
loglog(r_arr(12:end), -psfout,'o','color',get_color(im));
errorbar(r_arr(12:end), psfout, psfout_err, '.','color',get_color(im));
errorbar(r_arr(12:end), -psfout, psfout_err,'o','color',get_color(im));

ylim([3e-7,1.1])
xlim([4e-1,1.1e3])
grid on
xlabel('r [arcsec]', 'fontsize',15)
ylabel('PSF', 'fontsize',15)

psfdatallfields(ifield).im_in = im;
psfdatallfields(ifield).m_min_in = psfdat.m_min;
psfdatallfields(ifield).m_max_in = psfdat.m_max;
psfdatallfields(ifield).psf = [psfin, psfout];
psfdatallfields(ifield).psfps = [psfinps, psfoutps];
psfdatallfields(ifield).psf_err = [psfin_err, psfout_err];
psfdatallfields(ifield).psfps_err = [psfinps_err, psfoutps_err];
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
% save(sprintf('%s/psfdat',savedir),'psfdatallfields');
pltsavedir=(strcat(mypaths.alldat,'plots/TM',num2str(inst),'/'));
savename = sprintf('%s/%s_psf',pltsavedir,dt.name);
print(savename,'-dpng');close
end
end
%% plot PSF of all fields

for inst=[1,2]
figure
leg = {};
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
    load(sprintf('%s/psfdat',savedir),'psfdatallfields');
    r_arr = psfdatallfields(ifield).r_arr;
    psf = psfdatallfields(ifield).psfps;
    loglog(r_arr,psf,'linewidth',1.5,'DisplayName',dt.name);hold on
end
h=legend('show','Location','southwest');
set(h,'fontsize',12)
legend boxoff
ylim([3e-7,1.1])
xlim([4e-1,1.1e3])
grid on
xlabel('r [arcsec]', 'fontsize',15)
ylabel('PSF', 'fontsize',15)
bandname = ['I','H'];
title(sprintf('%s band',bandname(inst)),'fontsize',15);

pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
savename = sprintf('%s/psf',pltsavedir);
print(savename,'-dpng');close
end

