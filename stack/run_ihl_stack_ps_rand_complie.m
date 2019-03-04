% calculate the mean and error of 100 random position stacking

flight=40030;
inst=1;
ifield=5;
nsim = 100;
dt=get_dark_times(flight,inst,ifield);
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

ncountfile=strcat(savedir,'ciber_ps/',dt.name,'_stackcounts.txt');
T=readtable(ncountfile);
m_min_arr=T{:,2}';
m_max_arr=T{:,3}';
counts_arr=T{:,5}';
countg_arr=T{:,7}';
N_arr = unique([counts_arr, countg_arr]);
N_arr = N_arr(find(N_arr>0));

savename = strcat(savedir,'bk_ps/',dt.name);


M = csvread(strcat(savename,'_randprof','_1_1.csv'));
r_arr = M(:,1)';

for iN=1:numel(N_arr)
    
    N = N_arr(iN);
    profcb = zeros(size(r_arr));
    profcb_err = zeros(size(r_arr));
    profps = zeros(size(r_arr));
    profps_err = zeros(size(r_arr));
    
    for iter=1:nsim
        M = csvread(strcat(savename,'_randprof',...
            '_',num2str(N),'_',num2str(iter),'.csv'));

        profcb = profcb + M(:,2)';
        profps = profps + M(:,4)';
        profcb_err = profcb_err + M(:,2)'.^2;
        profps_err = profps_err + M(:,4)'.^2;
    end
    
    profcb_err = profcb_err - profcb;
    profps_err = profps_err - profps;
    
    profcb = profcb./nsim;
    profps = profps./nsim;
    profcb_err = profcb_err./nsim;
    profps_err = profps_err./nsim;
    
    profrand(iN).N = N;
    profrand(iN).r_arr = r_arr;
    profrand(iN).profcb = profcb;
    profrand(iN).profcb_err = profcb_err;
    profrand(iN).profps = profps;
    profrand(iN).profps_err = profps_err;
end

save(strcat(savedir,'bk_ps/',dt.name,'_profrand'),'profrand');
%% plot some background profile

savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

figure
setwinsize(gcf,1000,600)

iNplot_arr = [15,numel(profrand)];
for i=1:numel(iNplot_arr)
    
    iN = iNplot_arr(i);
    N = profrand(iN).N;
    r_arr = profrand(iN).r_arr;
    profcb = profrand(iN).profcb;
    profcb_err = profrand(iN).profcb_err;
    profps = profrand(iN).profps;
    profps_err = profrand(iN).profps_err;
    
    if i==1
        r_arr = r_arr.*0.95;
    elseif i==3
        r_arr = r_arr.*1.05;
    end

    subplot(2,2,1)
    errorbar(r_arr,profcb,profcb_err,'.-','color',get_color(i),...
        'DisplayName',strcat(num2str(N), ' CIBER stacking'));hold on   
    title('CIBER stacking');
    
    subplot(2,2,2)
        loglog(r_arr,profcb,'.-','color',get_color(i));hold on
        errorbar(r_arr,profcb,profcb_err,'.-','color',get_color(i));
        loglog(r_arr,-profcb,'o-','color',get_color(i));hold on
        errorbar(r_arr,-profcb,profcb_err,'o-','color',get_color(i));
    title('CIBER stacking');
    
    subplot(2,2,3)
    errorbar(r_arr,profps,profps_err,'.-','color',get_color(i),...
       'DisplayName',strcat(num2str(N), ' PanSTARRS stacking'));hold on   
    title('PanSTARRS stacking');
    
    subplot(2,2,4)
        loglog(r_arr,profps,'.-','color',get_color(i));hold on
        errorbar(r_arr,profps,profps_err,'.-','color',get_color(i));
        loglog(r_arr,-profps,'o-','color',get_color(i));hold on
        errorbar(r_arr,-profps,profps_err,'o-','color',get_color(i));
    title('PanSTARRS stacking');

end
subplot(2,2,1)
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
plot([4e-1,7e2],[0,0],'k--')

subplot(2,2,2)
xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')

subplot(2,2,3)
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
plot([4e-1,7e2],[0,0],'k--')

subplot(2,2,4)
xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')

suptitle(strcat(dt.name,' background stacking'));
savename=strcat(savedir,dt.name,'_background');
print(savename,'-dpng');%close
