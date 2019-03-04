% Sum up the flux from 1x1, 3x3, 5x5 pixels around the bright stars
% to get the cal factor ADU/fr to nW/m2/sr. The results written (by hand) in 
% get_cal_apf2nWpm2ps.

flight = 40030;
inst = 1;
mypaths=get_paths(flight);
savedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

for ifield = 4:8

dt=get_dark_times(flight,inst,ifield);
if ifield == 5
    cbmap =  - stackmapdat(ifield).map_last10;

else
    cbmap =  - stackmapdat(ifield).map;
end

strmask = maskdat.mask(ifield).strmask;
strnum = maskdat.mask(ifield).strnum;

mask_inst0 = stackmapdat(ifield).mask;
totmask = mask_inst0.*strmask;
sigmask1 = sigclip_mask(cbmap,totmask,3,5);
sm = fillpadsmooth(cbmap,sigmask1,2);
sigmask2 = sigclip_mask(sm,sigmask1,3,5);
sig_sp = find((totmask-sigmask2)==1);
mask_inst = mask_inst0;
mask_inst(sig_sp)=0;

m_binedges = 9:0.5:17;
m_min_arr = m_binedges(1:end-1);
m_max_arr = m_binedges(2:end);

I1_vec = zeros(size(m_min_arr));
I3_vec = zeros(size(m_min_arr));
I5_vec = zeros(size(m_min_arr));
IT_vec = zeros(size(m_min_arr));
I1_err = zeros(size(m_min_arr));
I3_err = zeros(size(m_min_arr));
I5_err = zeros(size(m_min_arr));
IT_err = zeros(size(m_min_arr));

for im = 1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    [I1, I3, I5, IT]=photometry_sum_2m...
        (flight,inst,ifield,m_min,m_max,cbmap,mask_inst,strmask,strnum);
    
    sp = find((I1 > median(I1) - 5 * std(I1)) & (I1 < median(I1) + 5 * std(I1))...
        & (I3 > median(I1) - 5 * std(I3)) & (I1 < median(I3) + 5 * std(I3)) ...
        & (I5 > median(I1) - 5 * std(I5)) & (I1 < median(I5) + 5 * std(I5)));
    
    I1_vec(im) = mean(I1(sp));
    I3_vec(im) = mean(I3(sp));
    I5_vec(im) = mean(I5(sp));
    IT_vec(im) = mean(IT(sp));
    I1_err(im) = std(I1(sp)) / sqrt(numel(sp));
    I3_err(im) = std(I3(sp)) / sqrt(numel(sp));
    I5_err(im) = std(I5(sp)) / sqrt(numel(sp));
    IT_err(im) = std(IT(sp)) / sqrt(numel(sp));   
end


if inst == 1
    prior = 300;
    fitmin = 1e5;
    fitmax = 1e6;
    if ifield == 5
        fitmax = 6e5;
    end
else
    prior = 100;
    fitmin = 3e4;
    fitmax = 2e5;
end

model =  @(p,x) p(1) .* x;

sp = find(IT_vec > fitmin & IT_vec < fitmax);
xData = I3_vec(sp);
yData = IT_vec(sp);
dxData = I3_err(sp);
dyData = IT_err(sp);

[fit.params,fit.dParams,fit.gof,fit.stddev] = ...
    fitChiSquare(xData,yData,model,prior,dxData,dyData);

figure
loglog(IT_vec, IT_vec./fit.params, 'k');hold on
errorbar(IT_vec, I1_vec, I1_err, I1_err, IT_err, IT_err, '.');
errorbar(IT_vec, I3_vec, I3_err, I3_err, IT_err, IT_err, '.');
errorbar(IT_vec, I5_vec, I5_err, I5_err, IT_err, IT_err, '.');
vline(fitmin);
vline(fitmax);
title(sprintf('TM%d, %s, cal = %.2f',inst, dt.name, fit.params));
xlabel('$\nu I_{\nu}[nW/m^2/sr]$', 'interpreter', 'latex');
ylabel('- ADU / fr');
drawnow
end