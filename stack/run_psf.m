%%%% run bright source stack for a range of mag bins %%%%%%
flight = 40030;
inst = 1;

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
m_min_arr = [4,4,4,4,12,13,15,16:19];
m_max_arr = [9,10,11,12,13,14,16,17:20];
%%
for ifield = 4:8
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
    
    if m_min<16
        Njack = 10;
        srcdat = tm_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
        'sample_type','jack_random','Nsub',Njack);

        if srcdat.Ns==0
            continue
        end

        [clipmaxs, clipmins, r_arr]=...
        stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
        mask_inst,strnum,1000,verbose,[],nan,true);
    else
        Njack = 16;
        srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
        'sample_type','jack_random','Nsub',Njack);
        % only stack 1000 src to speed up
        if srcdat.Ns>1000
            Nstot = 0;
            for isub=1:Njack
                Nisub = round(srcdat.sub(isub).Ns*1000/srcdat.Ns);
                srcdat.sub(isub).xs_arr = srcdat.sub(isub).xs_arr(1:Nisub);
                srcdat.sub(isub).ys_arr = srcdat.sub(isub).ys_arr(1:Nisub);
                srcdat.sub(isub).ms_arr = srcdat.sub(isub).ms_arr(1:Nisub);
                srcdat.sub(isub).Ns = Nisub;
                Nstot = Nstot + Nisub;
            end
            srcdat.Ns = Nstot;
        end
        [clipmaxs, clipmins, r_arr]=...
        stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
        mask_inst,strnum,1000,verbose,[],nan,false);
    end

    psfdat.r_arr = r_arr;
    mask_inst = squeeze(mask_inst(inst,:,:));
    
    for isub=1:Njack
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
        psfdat.sub(isub).profcbs = profcbs;
        psfdat.sub(isub).profpss = profpss;        
        psfdat.sub(isub).profhits = profhits;
    end
    
    %%% profile combining all subset
    profcbs = zeros(size(psfdat.r_arr));
    profpss = zeros(size(psfdat.r_arr));
    profhits = zeros(size(psfdat.r_arr));
    counts = 0;
    for isub=1:Njack
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
    for isub=1:Njack
        jackcbs = profcbs - psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        jackpss = profpss - psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        jackhits = profhits - psfdat.sub(isub).profhits;
        psfdat.jack(isub).profcbs = jackcbs./jackhits;
        psfdat.jack(isub).profpss = jackpss./jackhits;
    end
    
    %%% error bar with jackknife
    errcbs = zeros(size(psfdat.r_arr));
    errpss = zeros(size(psfdat.r_arr));
    for isub=1:Njack
        errcbs = errcbs + ...
            (psfdat.jack(isub).profcbs - psfdat.all.profcbs).^2;        
        errpss = errpss + ...
            (psfdat.jack(isub).profpss - psfdat.all.profpss).^2;
    end
    psfdat.errjack.profcbs = sqrt(errcbs.*((Njack-1)/Njack));
    psfdat.errjack.profpss = sqrt(errpss.*((Njack-1)/Njack));
    
    psfdatall(im).psfdat = psfdat;
end
psfdatallfields(ifield).psfdatall = psfdatall;
end
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/psfdat',savedir),'psfdatallfields');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
flight = 40030;
mypaths=get_paths(flight);
for inst=[1,2]
    savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
    load(sprintf('%s/psfdat',savedir),'psfdatallfields');
    psfdatallfieldsold = psfdatallfields;
    clear psfdatallfields;
    for ifield=4:8
        dt=get_dark_times(flight,inst,ifield);
        psfdatallfields(ifield).psfdatall = psfdatallfieldsold(ifield).psfdatall;
    end
    save(sprintf('%s/psfdat',savedir),'psfdatallfields');
end
%% plot the stack and get the combined PSF
flight = 40030;

mypaths=get_paths(flight);
m_min_arr = [4,4,4,4,12,13,15,16:19];
m_max_arr = [9,10,11,12,13,14,16,17:20];

for ifield=4:8
for inst=[1,2]
dt=get_dark_times(flight,inst,ifield);
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
vline(r_arr(9),'k--');
vline(r_arr(11),'k--');
h=legend('show','Location','northeast');
set(h,'fontsize',6)
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
for imin=8:11
plotoff = 3^-(imin-8);
psfdat = psfdatall(imin).psfdat;
norm = 1./psfdat.all.profcbs(1);
psfin = psfdat.all.profcbs(1:9).*norm;
psfin_err = psfdat.errjack.profcbs(1:9).*norm;
normps = 1./psfdat.all.profpss(1);
psfinps = psfdat.all.profpss(1:9).*normps;
psfinps_err = psfdat.errjack.profpss(1:9).*normps;

loglog(r_arr(1:9), psfin.*plotoff,'.-','color',get_color(imin));hold on
loglog(r_arr(1:9), -psfin.*plotoff,'o','color',get_color(imin));
errorbar(r_arr(1:9), psfin.*plotoff, psfin_err.*plotoff, ...
    '.','color',get_color(imin));
errorbar(r_arr(1:9), -psfin.*plotoff, psfin_err.*plotoff,...
    'o','color',get_color(imin));

immid = 6;
psfdat = psfdatall(immid).psfdat;
norm = psfin(9)./psfdat.all.profcbs(9);
psfmid = psfdat.all.profcbs(9:11).*norm;
psfmid_err = psfdat.errjack.profcbs(9:11).*norm;
normps = psfinps(9)./psfdat.all.profpss(9);
psfmidps = psfdat.all.profpss(9:11).*normps;
psfmidps_err = psfdat.errjack.profpss(9:11).*normps;

loglog(r_arr(9:11), psfmid.*plotoff,'.-','color',get_color(immid));hold on
loglog(r_arr(9:11), -psfmid.*plotoff,'o','color',get_color(immid));
errorbar(r_arr(9:11), psfmid.*plotoff, psfmid_err.*plotoff, ...
    '.','color',get_color(immid));
errorbar(r_arr(9:11), -psfmid.*plotoff, psfmid_err.*plotoff,...
    'o','color',get_color(immid));

[~,imout] = max(sum(snrs'));
psfdat = psfdatall(imout).psfdat;
norm = psfmid(end)./psfdat.all.profcbs(11);
psfout = psfdat.all.profcbs(11:end).*norm;
psfout_err = psfdat.errjack.profcbs(11:end).*norm;
normps = psfmidps(end)./psfdat.all.profpss(11);
psfoutps = psfdat.all.profpss(11:end).*normps;
psfoutps_err = psfdat.errjack.profpss(11:end).*normps;

psfdatallfields(ifield).psf(imin-7).r_arr = r_arr;
psfdatallfields(ifield).psf(imin-7).idx_r_cut = [9,11];
psfdatallfields(ifield).psf(imin-7).r_cut = r_arr([9,11]);
psfdatallfields(ifield).psf(imin-7).im_in = imin;
psfdatallfields(ifield).psf(imin-7).m_min_in = psfdatall(imin).psfdat.m_min;
psfdatallfields(ifield).psf(imin-7).m_max_in = psfdatall(imin).psfdat.m_max;
psfdatallfields(ifield).psf(imin-7).im_mid = immid;
psfdatallfields(ifield).psf(imin-7).m_min_mid = psfdatall(immid).psfdat.m_min;
psfdatallfields(ifield).psf(imin-7).m_max_mid = psfdatall(immid).psfdat.m_max;
psfdatallfields(ifield).psf(imin-7).im_out = imout;
psfdatallfields(ifield).psf(imin-7).m_min_out = psfdatall(imout).psfdat.m_min;
psfdatallfields(ifield).psf(imin-7).m_max_out = psfdatall(imout).psfdat.m_max;
psfdatallfields(ifield).psf(imin-7).psf = [psfin, psfmid(2:end), psfout(2:end)];
psfdatallfields(ifield).psf(imin-7).psfps = ...
    [psfinps, psfmidps(2:end), psfoutps(2:end)];
psfdatallfields(ifield).psf(imin-7).psf_err = ...
    [psfin_err, psfmid_err(2:end), psfout_err(2:end)];
psfdatallfields(ifield).psf(imin-7).psfps_err = ...
    [psfinps_err, psfmidps_err(2:end), psfoutps_err(2:end)];

loglog(r_arr(11:end), psfout.*plotoff,'.-','color',get_color(imout));hold on
loglog(r_arr(11:end), -psfout.*plotoff,'o','color',get_color(imout));
errorbar(r_arr(11:end), psfout.*plotoff, psfout_err.*plotoff, ...
    '.','color',get_color(imout));
errorbar(r_arr(11:end), -psfout.*plotoff, psfout_err.*plotoff,...
    'o','color',get_color(imout));
end
ylim([3e-7,1.1])
xlim([4e-1,1.1e3])
grid on
xlabel('r [arcsec]', 'fontsize',15)
ylabel('PSF', 'fontsize',15)

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/psfdat',savedir),'psfdatallfields');
pltsavedir=(strcat(mypaths.alldat,'plots/TM',num2str(inst),'/'));
savename = sprintf('%s/%s_psf',pltsavedir,dt.name);
print(savename,'-dpng');close
end
end
%% plot PSF of all fields
flight = 40030;
mypaths=get_paths(flight);

for inst=[1,2]
figure
leg = {};
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
    load(sprintf('%s/psfdat',savedir),'psfdatallfields');
    r_arr = psfdatallfields(ifield).psf(1).r_arr;
    psf = psfdatallfields(ifield).psf(1).psf;
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

pltsavedir=(strcat(mypaths.alldat,'plots/TM',num2str(inst),'/'));
savename = sprintf('%s/psf',pltsavedir);
% print(savename,'-dpng');close
end
