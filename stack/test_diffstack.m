flight=40030;
inst=1;
sample_type = 'jack_random';
subpix=false;
rmin=nan;
mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end
m_min_arr = stackmapdat(4).m_min_arr;
m_max_arr = stackmapdat(4).m_max_arr;
m_min_arr = m_min_arr(10:13);
m_max_arr = m_max_arr(10:13);

dx = 1200;
verbose = 0;

ifield1=4;
ifield2=5;
%% get FF for field1 - field2
FF = zeros(1024);
FFhit = zeros(1024);
for ifield=4:8
    if ifield~=ifield1 & ifield~=ifield2
    map = stackmapdat(ifield).rawmap.*cal(ifield).apf2nWpm2ps;
    strmask = stackmapdat(ifield).strmask;
    mask_inst = stackmapdat(ifield).mask;
    Q1 = quantile(map(find(mask_inst.*strmask)),0.25);
    Q3 = quantile(map(find(mask_inst.*strmask)),0.75);
    IQR = Q3-Q1;
    clipmin = Q1 - 3*IQR;
    clipmax = Q3 + 3*IQR;
    sigmask = mask_inst.*strmask;
    sigmask((map>clipmax) | (map<clipmin)) = 0;
    
    FF = FF + map.*sigmask./sqrt(mean(map(find(sigmask))));
    FFhit = FFhit + sigmask.*sqrt(mean(map(find(sigmask))));
    end
end
FF = FF./FFhit; FF(find(FFhit==0))=0;
%% process the field diff map
cal = get_cal_apf2nWpm2ps(inst);
rawmap5 = stackmapdat(ifield1).rawmap.*cal(ifield1).apf2nWpm2ps;
rawmap6 = stackmapdat(ifield2).rawmap.*cal(ifield2).apf2nWpm2ps;
dmap = rawmap5 - rawmap6;
dstrmask = stackmapdat(ifield1).strmask.*stackmapdat(ifield2).strmask;
mask_inst = stackmapdat(4).mask;
Q1 = quantile(dmap(find(mask_inst.*dstrmask)),0.25);
Q3 = quantile(dmap(find(mask_inst.*dstrmask)),0.75);
IQR = Q3-Q1;
clipmin = Q1 - 3*IQR;
clipmax = Q3 + 3*IQR;
sigmask = mask_inst.*dstrmask;
sigmask((dmap>clipmax) | (dmap<clipmin) | (FF==0)) = 0;
sig_sp = find((mask_inst.*dstrmask-sigmask)==1);
mask_inst_clip = mask_inst;
mask_inst_clip(sig_sp)=0;
dmap = dmap./FF;
dmap(find(FF==0))=0;
dmap = dmap - mean(dmap(find(sigmask)));
%%
dcbmap = dmap;
dpsmap = stackmapdat(ifield1).psmap - stackmapdat(ifield2).psmap;
dstrmask = stackmapdat(ifield1).strmask.*stackmapdat(ifield2).strmask;
dstrnum = stackmapdat(ifield1).strnum + stackmapdat(ifield2).strnum;
dmask_inst = mask_inst_clip;
%%
figure
setwinsize(gcf,1200,300)
subplot(1,3,1)
imageclip(stackmapdat(ifield1).cbmap.*stackmapdat(ifield1).strmask.*mask_inst);
title('elat10');
subplot(1,3,2)
imageclip(stackmapdat(ifield2).cbmap.*stackmapdat(ifield2).strmask.*mask_inst);
title('elat30');
subplot(1,3,3)
imageclip(dcbmap.*dstrmask.*dmask_inst);
title('elat10 - elat30');
%% 
ifield=ifield1;
dt=get_dark_times(flight,inst,ifield);

cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;

psmap = dcbmap;
mask_inst = mask_inst.*dmask_inst;
strmask = dstrmask;
strnum = dstrnum;
%%
for im= 3:4

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    stackdat.m_min = m_min;
    stackdat.m_max = m_max;
    
    srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,stackband,...
    'sample_type',sample_type);

    [clipmaxs, clipmins, rbins]=...
    stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
    mask_inst,strnum,1000,verbose,stackband,[],rmin);
    stackdat.r_arr = rbins;

    for isub=1:16
        [stampercbs,stamperpss,hitmaps,profcbs,profpss,profhits] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,rmin,clipmaxs,clipmins,...
            srcdat.sub(isub).xs_arr,srcdat.sub(isub).ys_arr,...
            srcdat.sub(isub).ms_arr,subpix);

        [stampercbg,stamperpsg,hitmapg,profcbg,profpsg,profhitg] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,rmin,clipmaxs,clipmins,...
            srcdat.sub(isub).xg_arr,srcdat.sub(isub).yg_arr,...
            srcdat.sub(isub).mg_arr,subpix);

        fprintf('stack %s, %d<m<%d, %s, isub %d\n',...
            dt.name,m_min,m_max,sample_type,isub);

        stackdat.sub(isub).counts = srcdat.sub(isub).Ns;
        stackdat.sub(isub).countg = srcdat.sub(isub).Ng;
        stackdat.sub(isub).profcbs = profcbs;
        stackdat.sub(isub).profcbg = profcbg;
        stackdat.sub(isub).profpss = profpss;
        stackdat.sub(isub).profpsg = profpsg;
        stackdat.sub(isub).profhits = profhits;
        stackdat.sub(isub).profhitg = profhitg;
    end
   
    sp100 = find(stackdat.r_arr>100);
  %% profile combining all subset
    profcbs = zeros(size(stackdat.r_arr));
    profcbg = zeros(size(stackdat.r_arr));
    profpss = zeros(size(stackdat.r_arr));
    profpsg = zeros(size(stackdat.r_arr));
    profhits = zeros(size(stackdat.r_arr));
    profhitg = zeros(size(stackdat.r_arr));
    counts = 0;
    countg = 0;
    for isub=1:16
        profcbs = profcbs + stackdat.sub(isub).profcbs.*stackdat.sub(isub).profhits;
        profcbg = profcbg + stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        profpss = profpss + stackdat.sub(isub).profpss.*stackdat.sub(isub).profhits;
        profpsg = profpsg + stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        profhits = profhits + stackdat.sub(isub).profhits;
        profhitg = profhitg + stackdat.sub(isub).profhitg;
        counts = counts + stackdat.sub(isub).counts;
        countg = countg + stackdat.sub(isub).countg;
    end
    stackdat.all.profcbs = profcbs./profhits;
    stackdat.all.profcbg = profcbg./profhitg;
    stackdat.all.profpss = profpss./profhits;
    stackdat.all.profpsg = profpsg./profhitg;
    stackdat.all.counts = counts;
    stackdat.all.countg = countg;
    
    stackdat.all.profcbs100 = sum(profcbs(sp100))./sum(profhits(sp100));
    stackdat.all.profcbg100 = sum(profcbg(sp100))./sum(profhitg(sp100));
    stackdat.all.profpss100 = sum(profpss(sp100))./sum(profhits(sp100));
    stackdat.all.profpsg100 = sum(profpsg(sp100))./sum(profhitg(sp100));

  %% profile of jackknife samples (leave one out)
    for isub=1:16
        jackcbs = profcbs - stackdat.sub(isub).profcbs.*stackdat.sub(isub).profhits;
        jackcbg = profcbg - stackdat.sub(isub).profcbg.*stackdat.sub(isub).profhitg;
        jackpss = profpss - stackdat.sub(isub).profpss.*stackdat.sub(isub).profhits;
        jackpsg = profpsg - stackdat.sub(isub).profpsg.*stackdat.sub(isub).profhitg;
        jackhits = profhits - stackdat.sub(isub).profhits;
        jackhitg = profhitg - stackdat.sub(isub).profhitg;
        stackdat.jack(isub).profcbs = jackcbs./jackhits;
        stackdat.jack(isub).profcbg = jackcbg./jackhitg;
        stackdat.jack(isub).profpss = jackpss./jackhits;
        stackdat.jack(isub).profpsg = jackpsg./jackhitg; 

        stackdat.jack(isub).profcbs100 = sum(jackcbs(sp100))./sum(jackhits);
        stackdat.jack(isub).profcbg100 = sum(jackcbg(sp100))./sum(jackhitg);
        stackdat.jack(isub).profpss100 = sum(jackpss(sp100))./sum(jackhits);
        stackdat.jack(isub).profpsg100 = sum(jackpsg(sp100))./sum(jackhitg); 
        
    end
  %% error bar with jackknife
    errcbs = zeros(size(stackdat.r_arr));
    errcbg = zeros(size(stackdat.r_arr));
    errpss = zeros(size(stackdat.r_arr));
    errpsg = zeros(size(stackdat.r_arr));

    errcbs100 = 0;
    errcbg100 = 0;
    errpss100 = 0;
    errpsg100 = 0;

    for isub=1:16
    errcbs = errcbs + (stackdat.jack(isub).profcbs - stackdat.all.profcbs).^2;
    errcbg = errcbg + (stackdat.jack(isub).profcbg - stackdat.all.profcbg).^2;
    errpss = errpss + (stackdat.jack(isub).profpss - stackdat.all.profpss).^2;
    errpsg = errpsg + (stackdat.jack(isub).profpsg - stackdat.all.profpsg).^2;

    errcbs100 = errcbs100 + ...
        (stackdat.jack(isub).profcbs100 - stackdat.all.profcbs100).^2;
    errcbg100 = errcbg100 + ...
        (stackdat.jack(isub).profcbg100 - stackdat.all.profcbg100).^2;
    errpss100 = errpss100 + ...
        (stackdat.jack(isub).profpss100 - stackdat.all.profpss100).^2;
    errpsg100 = errpsg100 + ...
        (stackdat.jack(isub).profpsg100 - stackdat.all.profpsg100).^2;
    end
    stackdat.errjack.profcbs = sqrt(errcbs.*(15/16));
    stackdat.errjack.profcbg = sqrt(errcbg.*(15/16));
    stackdat.errjack.profpss = sqrt(errpss.*(15/16));
    stackdat.errjack.profpsg = sqrt(errpsg.*(15/16));
    stackdat.errjack.profcbs100 = sqrt(errcbs100.*(15/16));
    stackdat.errjack.profcbg100 = sqrt(errcbg100.*(15/16));
    stackdat.errjack.profpss100 = sqrt(errpss100.*(15/16));
    stackdat.errjack.profpsg100 = sqrt(errpsg100.*(15/16));
%%
    dat(im).stackdat = stackdat;
end
%%
r_arr = stackdat.r_arr;
figure
setwinsize(gcf,800,600)
for im=3:4
    m_min = dat(im).stackdat.m_min;
    m_max = dat(im).stackdat.m_max;

    profcbs = dat(im).stackdat.all.profcbs;
    profcbg = dat(im).stackdat.all.profcbg;
    errcbs = dat(im).stackdat.errjack.profcbs;
    errcbg = dat(im).stackdat.errjack.profcbg;

    profpss = dat(im).stackdat.all.profpss;
    profpsg = dat(im).stackdat.all.profpsg;
    errpss = dat(im).stackdat.errjack.profpss;
    errpsg = dat(im).stackdat.errjack.profpsg;

    subplot(2,2,im-2)
    semilogx(r_arr.*1.01, profcbs, 'r');hold on
    semilogx(r_arr.*0.99, profcbg, 'b');
    legend({'stars','galaxies'});
    hline(0,'k--');
    errorbar(r_arr.*1.01, profcbs, errcbs,'r.');
    errorbar(r_arr.*0.99, profcbg, errcbg,'b.');
    title(sprintf('%d<m<%d elat10',m_min,m_max));
    xlim([1e1,1.1e3])
    ylim([-5,6])

    subplot(2,2,im)
    semilogx(r_arr.*1.01, profpss, 'r');hold on
    semilogx(r_arr.*0.99, profpsg, 'b');
    legend({'stars','galaxies'});
    hline(0,'k--');
    errorbar(r_arr.*1.01, profpss, errpss,'r.');
    errorbar(r_arr.*0.99, profpsg, errpsg,'b.');
    title(sprintf('%d<m<%d elat10 - elat30',m_min,m_max));
    xlim([1e1,1.1e3])
    ylim([-5,6])
end
    print(sprintf('/Users/ytcheng/Desktop/ciber/fluc_test'),'-dpng');
