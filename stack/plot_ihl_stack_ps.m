%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ihl stacking maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
ifield=4;
npix=1200;

dt=get_dark_times(flight,inst,ifield);
loaddir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/');
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

load(strcat(loaddir,'bk_ps/',dt.name,'_profrand'));

m_min_arr = [0,8:22];
m_max_arr = [8:23];
%% get the PSF
psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

bestparam = fitpsfdat(ifield).bestparam;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);

radmap = make_radius_map(zeros(2*npix+1),npix,npix).*0.7;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);

profile = radial_prof(psfmap,ones(2*npix+1),npix+1,npix+1,1,25);
r_arr=profile.r*0.7;
profpsf_arr=(profile.prof)./profile.prof(1);
%% get the number counts data
ncountfile=strcat(loaddir,'ciber_ps/',dt.name,'_stackcounts.txt');
T=readtable(ncountfile);
mbot_arr=T{:,2}';
counts_arr=T{:,5}';
countg_arr=T{:,7}';
countstot_arr=T{:,4}';
countgtot_arr=T{:,6}';
%%
mask_inst_clip = ...
    fits_read(strcat(loaddir,'ciber_ps/',dt.name,'_mask_inst_clip.fits'));
frac = sum(mask_inst_clip(:))/1024/1024;

figure
semilogy(mbot_arr+0.5, countstot_arr./4./frac,'ro-');hold on
semilogy(mbot_arr+0.5, countgtot_arr./4./frac,'bo-');
semilogy(mbot_arr+0.5, counts_arr./4./frac,'ro--');
semilogy(mbot_arr+0.5, countg_arr./4./frac,'bo--');

ylim([1e0,1e5])
xlabel('m_{AB} (y band)')
ylabel('dN/dm_{AB}/deg^2')
legend({'stars','galaxies','stars stacked','galaxies stacked'});
title(dt.name)
savename=strcat(savedir,dt.name,'_hist');
print(savename,'-dpng');%close
%% stacking source with non-FF corrected map
%{
loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
num2str(inst),'/'));
m=19;
counts=sum(counts_arr(find(mbot_arr==m)));
countg=sum(countg_arr(find(mbot_arr==m)));
for itype=1:2
    if itype==1
        type='ciber_ps';        
    else
        type='panstarrs_ps';
    end
    stacks=fitsread(strcat(loaddir,'/',char(type),'/',dt.name,...
        '_stampers',num2str(m),'1.fits'));
    stackg=fitsread(strcat(loaddir,'/',char(type),'/',dt.name,...
        '_stamperg',num2str(m),'1.fits'));                    
    mstacks=fitsread(strcat(loaddir,'/ciber_ps/',dt.name,...
        '_hitmaps',num2str(m),'1.fits'));
    mstackg=fitsread(strcat(loaddir,'/ciber_ps/',dt.name,...
        '_hitmapg',num2str(m),'1.fits'));
    
    if itype==1
        stackscb=stacks./mstacks;
        stackgcb=stackg./mstackg;
    else
        stacks2m=stacks./mstacks;
        stackg2m=stackg./mstackg;
    end
    
end
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(stackscb(npix-1000:npix+1000,npix-1000:npix+1000));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,2)
imageclip(stackgcb(npix-1000:npix+1000,npix-1000:npix+1000));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
subplot(2,2,3)
imageclip(stacks2m(npix-1000:npix+1000,npix-1000:npix+1000));
title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,4)
imageclip(stackg2m(npix-1000:npix+1000,npix-1000:npix+1000));
title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

suptitle(strcat('PanSTARRS',{' '},num2str(m),'<mAB(y band)<',num2str(m+1)));

savename=strcat(savedir,dt.name,'_stackmaps_allquad',num2str(m),'largeFF');
print(savename,'-dpng');close
%}
%% plot the stacking map
for im = 10:13%1:numel(m_min_arr)
m_min = m_min_arr(im);
m_max = m_max_arr(im);
counts=sum(counts_arr(find(mbot_arr==m_min)));
countg=sum(countg_arr(find(mbot_arr==m_min)));

for itype=1:2
    if itype==1
        type='ciber_ps';        
    else
        type='panstarrs_ps';
    end
    stacks=fitsread(strcat(loaddir,'/',char(type),'/',dt.name,...
        '_stampers',num2str(m_min),'_',num2str(m_max),'.fits'));
    stackg=fitsread(strcat(loaddir,'/',char(type),'/',dt.name,...
        '_stamperg',num2str(m_min),'_',num2str(m_max),'.fits'));                    
    mstacks=fitsread(strcat(loaddir,'/ciber_ps/',dt.name,...
        '_hitmaps',num2str(m_min),'_',num2str(m_max),'.fits'));
    mstackg=fitsread(strcat(loaddir,'/ciber_ps/',dt.name,...
        '_hitmapg',num2str(m_min),'_',num2str(m_max),'.fits'));

    if itype==1
        stackscb=stacks./mstacks;
        stackgcb=stackg./mstackg;
    else
        stacksps=stacks./mstacks;
        stackgps=stackg./mstackg;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot stacking maps
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(stackscb(npix-100:npix+100,npix-100:npix+100));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,2)
imageclip(stackgcb(npix-100:npix+100,npix-100:npix+100));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
subplot(2,2,3)
imageclip(stacksps(npix-100:npix+100,npix-100:npix+100));
title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,4)
imageclip(stackgps(npix-100:npix+100,npix-100:npix+100));
title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

suptitle(strcat('PanSTARRS',{' '},num2str(m_min),'<mAB(y band)<',num2str(m_max)));

savename=strcat(savedir,dt.name,'_stackmaps',...
    num2str(m_min),'_',num2str(m_max));
print(savename,'-dpng');close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get the background stacking systematic profile

if counts~=0
    idx = find([profrand.N]==counts);
    profcb_bks = profrand(idx).profcb;
    profcb_err_bks = profrand(idx).profcb_err;
    profps_bks = profrand(idx).profps;
    profps_err_bks = profrand(idx).profps_err;
else
    profcb_bks = 0;
    profcb_err_bks = 0;
    profps_bks = 0;
    profps_err_bks = 0;    
end

if countg~=0
    idx = find([profrand.N]==countg);
    profcb_bkg = profrand(idx).profcb;
    profcb_err_bkg = profrand(idx).profcb_err;
    profps_bkg = profrand(idx).profps;
    profps_err_bkg = profrand(idx).profps_err;
else
    profcb_bkg = 0;
    profcb_err_bkg = 0;
    profps_bkg = 0;
    profps_err_bkg = 0;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot radial profile (unnormalizaed)

figure
setwinsize(gcf,1000,800)

ymin=[];
for itype=1:2
    if itype==1
        stacks=stackscb;
        stackg=stackgcb;
        prof_bks = profcb_bks;
        prof_err_bks = profcb_err_bks;
        prof_bkg = profcb_bkg;
        prof_err_bkg = profcb_err_bkg;       
    else
        stacks=stacksps;
        stackg=stackgps;
        prof_bks = profps_bks;
        prof_err_bks = profps_err_bks;
        prof_bkg = profps_bkg;
        prof_err_bkg = profps_err_bkg;
    end
    
    profile = radial_prof(stacks,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',3);
    r_arr=profile.r*0.7;
    profs_arr=profile.prof - prof_bks;
    errs_arr=sqrt(profile.err.^2 + prof_err_bks.^2);

    profile = radial_prof(stackg,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',3);
    profg_arr=profile.prof - prof_bkg;
    errg_arr=sqrt(profile.err.^2 + prof_err_bkg.^2);
    
    ymin=[ymin profs_arr profg_arr];
        
    subplot(2,2,1)
    if itype==1
        errorbar(r_arr,profs_arr,errs_arr,'r.-','DisplayName','CBstars');hold on
        errorbar(r_arr,profg_arr,errg_arr,'b.-','DisplayName','CBgals');
    else
        errorbar(r_arr,profs_arr,errs_arr,'m.-','DisplayName','PSstars');
        errorbar(r_arr,profg_arr,errg_arr,'c.-','DisplayName','PSgals');     
    end
    
    subplot(2,2,2)
    if itype==1
        loglog(r_arr.*0.98,profs_arr,'r.','markersize',10);hold on
        loglog(r_arr.*1.02,profg_arr,'b.','markersize',10);
        errorbar(r_arr.*0.98,profs_arr,errs_arr,'r.','markersize',10);
        errorbar(r_arr.*1.02,profg_arr,errg_arr,'b.','markersize',10);
        loglog(r_arr.*0.98,-profs_arr,'ro','markersize',5);hold on
        loglog(r_arr.*1.02,-profg_arr,'bo','markersize',5);
        errorbar(r_arr.*0.98,-profs_arr,errs_arr,'ro','markersize',5);
        errorbar(r_arr.*1.02,-profg_arr,errg_arr,'bo','markersize',5);
    else
        loglog(r_arr.*0.98,profs_arr,'m.','markersize',10);hold on
        loglog(r_arr.*1.02,profg_arr,'c.','markersize',10);
        errorbar(r_arr.*0.98,profs_arr,errs_arr,'m.','markersize',10);
        errorbar(r_arr.*1.02,profg_arr,errg_arr,'c.','markersize',10);        
        loglog(r_arr.*0.98,-profs_arr,'mo','markersize',5);hold on
        loglog(r_arr.*1.02,-profg_arr,'co','markersize',5);
        errorbar(r_arr.*0.98,-profs_arr,errs_arr,'mo','markersize',5);
        errorbar(r_arr.*1.02,-profg_arr,errg_arr,'co','markersize',5);
    end

end

subplot(2,2,1)
xlabel('arcsec')
ylabel('<I_{stack}>')
title(strcat('$\frac{\sum_i \rm{stamp}}{\sum_i 1}$'...
    ,'(unnormalized, linear)'),...
    'interpreter','latex','fontsize',15)
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');

set(h,'fontsize',10)
legend boxoff

ymin=ymin(find(ymin>0));
ymin=floor(log10(min(ymin)));
ymin=10^ymin;
subplot(2,2,2)
xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')
title(strcat('$\frac{\sum_i \rm{stamp}}{\sum_i 1}$'...
    ,'(unnormalized, log)'),...
    'interpreter','latex','fontsize',15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot radial profile normalized

profstack(im).m_min = m_min;
profstack(im).m_max = m_max;
profstack(im).r_arr = r_arr;
profstack(im).Ns = counts;
profstack(im).Ng = countg;

ymin=[];
for itype=1:2
    if itype==1
        stacks=stackscb;
        stackg=stackgcb;
        prof_bks = profcb_bks;
        prof_err_bks = profcb_err_bks;
        prof_bkg = profcb_bkg;
        prof_err_bkg = profcb_err_bkg;       
    else
        stacks=stacksps;
        stackg=stackgps;
        prof_bks = profps_bks;
        prof_err_bks = profps_err_bks;
        prof_bkg = profps_bkg;
        prof_err_bkg = profps_err_bkg;
    end
    
    norm = stacks(npix+1,npix+1);  
    profile = radial_prof(stacks,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',3);
    r_arr=profile.r*0.7;
    profs_arr=(profile.prof - prof_bks)./norm;
    errs_arr=sqrt(profile.err.^2 + prof_err_bks.^2)./norm;

    norm = stackg(npix+1,npix+1);  
    profile = radial_prof(stackg,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',3);
    profg_arr=(profile.prof - prof_bks)./norm;
    errg_arr=sqrt(profile.err.^2 + prof_err_bkg.^2)./norm;
    
    if itype==1
        profstack(im).profscb = profs_arr;
        profstack(im).profscb_err = errs_arr;
        profstack(im).profgcb = profg_arr;
        profstack(im).profgcb_err = errg_arr;        
    else
        profstack(im).profsps = profs_arr;
        profstack(im).profsps_err = errs_arr;
        profstack(im).profgps = profg_arr;
        profstack(im).profgps_err = errg_arr;                
    end
    
    ymin=[ymin profs_arr profg_arr];
    subplot(2,2,3)
    if itype==1
        errorbar(r_arr,profs_arr,errs_arr,'r.-','DisplayName','CBstars');hold on
        errorbar(r_arr,profg_arr,errg_arr,'b.-','DisplayName','CBgals');
        else
        errorbar(r_arr,profs_arr,errs_arr,'m.-','DisplayName','PSstars');
        errorbar(r_arr,profg_arr,errg_arr,'c.-','DisplayName','PSgals');     
    end
    
    subplot(2,2,4)
    if itype==1
        loglog(r_arr.*0.98,profs_arr,'r.','markersize',10);hold on
        loglog(r_arr.*1.02,profg_arr,'b.','markersize',10);
        errorbar(r_arr.*0.98,profs_arr,errs_arr,'r.','markersize',10);
        errorbar(r_arr.*1.02,profg_arr,errg_arr,'b.','markersize',10);
        loglog(r_arr.*0.98,-profs_arr,'ro','markersize',5);hold on
        loglog(r_arr.*1.02,-profg_arr,'bo','markersize',5);
        errorbar(r_arr.*0.98,-profs_arr,errs_arr,'ro','markersize',5);
        errorbar(r_arr.*1.02,-profg_arr,errg_arr,'bo','markersize',5);
    else
        loglog(r_arr.*0.98,profs_arr,'m.','markersize',10);hold on
        loglog(r_arr.*1.02,profg_arr,'c.','markersize',10);
        errorbar(r_arr.*0.98,profs_arr,errs_arr,'m.','markersize',10);
        errorbar(r_arr.*1.02,profg_arr,errg_arr,'c.','markersize',10);        
        loglog(r_arr.*0.98,-profs_arr,'mo','markersize',5);hold on
        loglog(r_arr.*1.02,-profg_arr,'co','markersize',5);
        errorbar(r_arr.*0.98,-profs_arr,errs_arr,'mo','markersize',5);
        errorbar(r_arr.*1.02,-profg_arr,errg_arr,'co','markersize',5);
    end

end

subplot(2,2,3)
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
ylim([-0.2,1.4])
title(strcat('$\frac{\sum_i \rm{stamp}}{\sum_i 1}$'...
    ,'(normalized, linear)'),...
    'interpreter','latex','fontsize',15)
set(gca, 'XScale', 'log')
h=legend('show','Location','southwest');
text(100,1.3,strcat('stack',{' '}, num2str(counts),{' '},'stars'));
text(100,1.2,strcat('stack',{' '}, num2str(countg),{' '},'gals'));

set(h,'fontsize',10)
legend boxoff

ymin=ymin(find(ymin>0));
ymin=floor(log10(min(ymin)));
ymin=10^ymin;
subplot(2,2,4)
loglog(r_arr,profpsf_arr,'k--');
xlim([4e-1,1e3])
ylim([ymin,2e0])
title(strcat('$\frac{\sum_i \rm{stamp}}{\sum_i 1}$'...
    ,'(normalized, log)'),'interpreter','latex','fontsize',15)
xlabel('arcsec')
ylabel('<I_{stack}>')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
suptitle(strcat('PanSTARRS',{' '},num2str(m_min),'<mAB(y band)<',num2str(m_max)));
savename=strcat(savedir,dt.name,'_rprof',num2str(m_min),'_',num2str(m_max),'bksub');
print(savename,'-dpng');%close

end
save(sprintf('%s/ciber_ps/%s_profstack',loaddir,dt.name),'profstack');
%% plot unnorm profile of different mag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
figure
for itype=1:4
    if itype==1
        data_type = 'ciber';src_type = 'g';name = 'CIBER galaxies';c = 'b';
    elseif itype==2
        data_type = 'ciber';src_type = 's';name = 'CIBER stars';c = 'r';
    elseif itype==3
        data_type = 'panstarrs';src_type = 'g';name = 'Pan-STARRS galaxies';c = 'c';
    else
        data_type = 'panstarrs';src_type = 's';name = 'Pan-STARRS stars';c = 'm';
    end
        
for m=12:20
    loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

    stack=fitsread(strcat(loaddir,'/',char(data_type),'_ps/',dt.name,...
        '_stamper',src_type,num2str(m),'.fits'));                    
    mstack=fitsread(strcat(loaddir,'/ciber_ps/',dt.name,...
        '_hitmap',src_type,num2str(m),'.fits'));

    stack = stack./mstack;

    profile = radial_prof(stack,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',3);
    r_arr=profile.r*0.7;
    prof_arr=(profile.prof);
    err_arr=profile.err;
    
    if m==20
        plot(r_arr,prof_arr,'.-','color',c);hold on
    else
        plot(r_arr,prof_arr,'.-','color',c);hold on
    end
    drawnow
end
end

xlabel('arcsec')
ylabel('<I_{stack}>')
title(strcat('$\frac{\sum_i \rm{stamp}}{\sum_i 1}$'...
    ,'(unnormalized)'),...
    'interpreter','latex','fontsize',15)
%xlim([4e-1,7e2])
xlim([4e-1,1e3])
%ylim([1e-1,1e5])
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
%h=legend('show','Location','southwest');
%set(h,'fontsize',10)
%legend boxoff

savename=strcat(savedir,dt.name,'_rprof_allmag');
print(savename,'-dpng');%close
%}