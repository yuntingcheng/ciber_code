function ihl_stack_compile_sides(flight,inst,ifield)

mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

stackdir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/sides/'));

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat_sides',loaddir),'stackmapdat');

npix = 1200;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;

%%% get the stacking maps %%%
datacb = zeros(2*npix+1, 2*npix+1);
dataps = zeros(2*npix+1, 2*npix+1);
mdata = zeros(2*npix+1, 2*npix+1);
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
   
    datacb = datacb + fitsread(strcat(stackdir,dt.name,'_stampercb',...
        num2str(m_min),'_',num2str(m_max),'.fits'));
    dataps = dataps + fitsread(strcat(stackdir,dt.name,'_stamperps',...
        num2str(m_min),'_',num2str(m_max),'.fits'));
    mdata = mdata + fitsread(strcat(stackdir,dt.name,'_hitmap',...
        num2str(m_min),'_',num2str(m_max),'.fits'));
    
end
stackcb=datacb./mdata;
stackps=dataps./mdata;

%%% plot the stacking maps %%%
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(stackps);
caxis([-1.5,1.5])
title('sim map m_{AB}^y < 20');
subplot(2,2,2)
imageclip(stackcb);
caxis([-1.5,1.5])
title('sim map m_{AB}^y < 28');
subplot(2,2,3)
imageclip(stackps(npix-100:npix+100,npix-100:npix+100));
caxis([-1,2])
title('stack 16 < m_{AB}^y < 19');
subplot(2,2,4)
imageclip(stackcb(npix-100:npix+100,npix-100:npix+100));
caxis([-1,2])
title('stack 16 < m_{AB}^y < 19');
suptitle('SIDES sim stacking 16 < m_{AB}^y < 19');

savename=strcat(pltsavedir,dt.name,'_stackmaps_sides');
print(savename,'-dpng');%close

%%% get the stacking profile %%%
prof = radial_prof(stackcb,ones(2*npix+1),npix+1,npix+1,1,25,...
    'sig',3,'iter_clip',3);
profcb = prof.prof;
errcb = prof.err;
norm = stackcb(npix+1,npix+1);  
profcbnorm = profcb./norm;
errcbnorm = errcb./norm;

prof = radial_prof(stackps,ones(2*npix+1),npix+1,npix+1,1,25,...
    'sig',3,'iter_clip',3);
profps = prof.prof;
errps = prof.err;
norm = stackps(npix+1,npix+1);  
profpsnorm = profps./norm;
errpsnorm = errps./norm;


diff = profcbnorm - profpsnorm;
differr = sqrt(errcbnorm.^2 + errpsnorm.^2);

r_arr = prof.r.*0.7;

datadir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/excessdat',datadir),'excessdat');

figure
setwinsize(gcf,1500,300)


subplot(1,3,1)
errorbar(r_arr,profcbnorm,errcbnorm,'r.-',...
    'DisplayName','sim map m_{AB}^y < 28');hold on
errorbar(r_arr,profpsnorm,errcbnorm,'b.-',...
    'DisplayName','sim map m_{AB}^y < 20');

xlabel('arcsec')
ylabel('<I_{stack}>')
title('normalized stacking profile (linear)','fontsize',15)
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

subplot(1,3,2)
loglog(r_arr.*0.98,profcbnorm,'r.','markersize',10);hold on
errorbar(r_arr.*0.98,profcbnorm,errcbnorm,'r.','markersize',10);
errorbar(r_arr.*0.98,-profcbnorm,errcbnorm,'ro','markersize',5);
loglog(r_arr.*1.02,profpsnorm,'b.','markersize',10);hold on
errorbar(r_arr.*1.02,profpsnorm,errpsnorm,'b.','markersize',10);
errorbar(r_arr.*1.02,-profpsnorm,errpsnorm,'bo','markersize',5);

xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')
title('normalized stacking profile (log)','fontsize',15)

subplot(1,3,3)
loglog(r_arr,diff,'k.','markersize',10);hold on
loglog(excessdat.r_arr,excessdat.all_avg.prof,...
    'r.','markersize',10);hold on
h=legend({'SIDES diff','CB - PS excess'},'Location','northeast');
set(h,'fontsize',10)
legend boxoff
errorbar(r_arr,diff,differr,'k.','markersize',10);
errorbar(r_arr,-diff,differr,'ko','markersize',5);
errorbar(excessdat.r_arr,excessdat.all_avg.prof,excessdat.all_avg.err,...
    'r.','markersize',10);
errorbar(excessdat.r_arr,-excessdat.all_avg.prof,excessdat.all_avg.err,...
'ro','markersize',5);

xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('normalized <I_{stack}>')
title('diff','fontsize',15)

savename=strcat(pltsavedir,dt.name,'_prof_sides');
print(savename,'-dpng');%close

return