flight = 40030;
inst =2;
mypaths=get_paths(flight);
verbose = 0;

m_min_arr = [0,8:22];
m_max_arr = [8:23];

Icorrdat.m_min_arr = m_min_arr;
Icorrdat.m_max_arr = m_max_arr;

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat0',loaddir),'stackmapdat');

for ifield = 4:8

dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
if ifield == 5
    cbmap = stackmapdat(ifield).cbmap_last10;
end
psmap = stackmapdat(ifield).psmap;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;

corrg_arr = zeros(size(m_min_arr));
corrs_arr = zeros(size(m_min_arr));
corru_arr = zeros(size(m_min_arr));
for im = 1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    type = 1;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_ps_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum,'y');
    if sum(isnan(cliplim_arrin(:))) == 0
    [stampercbc,stamperpsc,hitmapc,~]=...
        stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,100,verbose,0,0,'y');
    
    profcb = radial_prof(stampercbc./hitmapc,ones(2*50+1),50+1,50+1,1,15,...
        0,'weight',hitmapc);
    profps = radial_prof(stamperpsc./hitmapc,ones(2*50+1),50+1,50+1,1,15,...
        0,'weight',hitmapc);
    corrg_arr(im) = mean(profcb.prof(1:5)./profps.prof(1:5));
    end
    
    type = -1;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_ps_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum,'y');
    if sum(isnan(cliplim_arrin(:))) == 0
    [stampercbc,stamperpsc,hitmapc,~]=...
        stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,100,verbose,0,0,'y');
    
    profcb = radial_prof(stampercbc./hitmapc,ones(2*50+1),50+1,50+1,1,15,...
        0,'weight',hitmapc);
    profps = radial_prof(stamperpsc./hitmapc,ones(2*50+1),50+1,50+1,1,15,...
        0,'weight',hitmapc);
    corrs_arr(im) = mean(profcb.prof(1:5)./profps.prof(1:5));
    end

    type = 2;
    [cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_ps_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum,'y');
    if sum(isnan(cliplim_arrin(:))) == 0
    [stampercbc,stamperpsc,hitmapc,~]=...
        stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,100,verbose,0,0,'y');
    
    profcb = radial_prof(stampercbc./hitmapc,ones(2*50+1),50+1,50+1,1,15,...
        0,'weight',hitmapc);
    profps = radial_prof(stamperpsc./hitmapc,ones(2*50+1),50+1,50+1,1,15,...
        0,'weight',hitmapc);
    corru_arr(im) = mean(profcb.prof(1:5)./profps.prof(1:5));
    end  
end
Icorrdat.field(ifield).corrg_arr = corrg_arr;
Icorrdat.field(ifield).corrs_arr = corrs_arr;
Icorrdat.field(ifield).corru_arr = corru_arr;

end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/Icorrdat',savedir),'Icorrdat');

%%
flight = 40030;
inst =1;
mypaths=get_paths(flight);

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/Icorrdat',savedir),'Icorrdat');

figure
setwinsize(gcf,1000,300)
ms_arr = [];
corrs_arr = [];
mg_arr = [];
corrg_arr = [];
for ifield =4:8
dt=get_dark_times(flight,inst,ifield);
corrg = [];
mg = [];
corrs = [];
ms = [];
corru = [];
mu = [];
for im = 1:numel(m_min_arr)
    m_min = Icorrdat.m_min_arr(im);
    m_max = Icorrdat.m_max_arr(im);
    m = (m_min + m_max)./2;
        if Icorrdat.field(ifield).corrg_arr(im) ~= 0
            mg = [mg, m];
            corrg = [corrg, Icorrdat.field(ifield).corrg_arr(im)];
        end
        if Icorrdat.field(ifield).corrs_arr(im) ~= 0
            ms = [ms, m];
            corrs = [corrs, Icorrdat.field(ifield).corrs_arr(im)];
        end
        if Icorrdat.field(ifield).corru_arr(im) ~= 0
            mu = [mu, m];
            corru = [corru, Icorrdat.field(ifield).corru_arr(im)];
        end
end
if ifield ~= 5
    ms_arr = [ms_arr, ms(ms > 14 & ms < 19)];
    corrs_arr = [corrs_arr, corrs(ms > 14 & ms < 19)];
    mg_arr = [mg_arr, mg(mg > 16 & mg < 20)];
    corrg_arr = [corrg_arr, corrg(mg > 16 & mg < 20)];
end
subplot(1,3,1)
plot(ms, corrs, '.-', 'color', get_color(ifield-3), ...
    'markersize',15, 'DisplayName', strcat(dt.name) ); hold on
subplot(1,3,2)
plot(mg, corrg, '.-', 'color', get_color(ifield-3), ...
    'markersize',15, 'DisplayName', strcat(dt.name) ); hold on
subplot(1,3,3)
plot(mu, corru, '.-', 'color', get_color(ifield-3), ...
    'markersize',15, 'DisplayName', strcat(dt.name) ); hold on
end

ps = polyfit(ms_arr, corrs_arr,0);
pg = polyfit(mg_arr, corrg_arr,1);

subplot(1,3,1)
plot(ms_arr, ps * ones(size(ms_arr)), 'k', 'linewidth',2, 'DisplayName', 'fit');
subplot(1,3,2)
plot(mg_arr, mg_arr * pg(1) + pg(2), 'k', 'linewidth',2, 'DisplayName', 'fit');
subplot(1,3,3)
plot(ms_arr, ps * ones(size(ms_arr)), 'k', 'linewidth',2, 'DisplayName', 'fit');

subplot(1,3,1)
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
xlabel('y band mag')
ylabel('stack CIBER / stack Sim')
title('stars')
ylim([0,6])

subplot(1,3,2)
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
xlabel('y band mag')
ylabel('stack CIBER / stack Sim')
title('galaxies')
ylim([0,6])

subplot(1,3,3)
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
xlabel('y band mag')
ylabel('stack CIBER / stack Sim')
title('undefined')
ylim([0,6])

fprintf('band %d star Icorr = %.2f * Ilin \n',inst, ps);
fprintf('band %d gals Icorr = (%.2f * my + %.2f) * Ilin \n',inst, pg(1),pg(2));