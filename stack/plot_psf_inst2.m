%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (TM2) Take the PSF map from 2MASS stacking,
% plot the PSF map, Bl, angular average radial profile.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
loaddir='/Users/ytcheng/ciber/doc/20170617_Stacking/psf/';
loaddir1='/Users/ytcheng/ciber/data/iband_psf_quads/';
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/plots/TM2/');
quad_arr=['A','B','C','D'];
%% plot the psf
figure
setwinsize(gcf,1000,1250)
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'yt/inst2/j0_14/',dt.name,'_',quad,'.fits');
    psfytc=fitsread(file);
    
    subplot(5,4,(ifield-4)*4+iquad)
    imageclip(psfytc(300:500,300:500));
    title(strcat(dt.name,'\_',quad))
    
end
end
savename=strcat(savedir,'psfmap');
print(savename,'-dpng');close

%% plot bl
figure
setwinsize(gcf,1500,500)

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
for iquad=1:4
    
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'yt/inst2/j0_14/',dt.name,'_',quad,'.fits');
    psfytc=fitsread(file);
    
    psfytc(find(psfytc~=psfytc))=0;
    
    [bly,lfine]=get_angular_spec(psfytc,psfytc,0.7,'nbins',30);

    lfine=lfine(find(bly));
    bly=bly(find(bly));
    bly=bly./mean(bly(1:3));
    
    subplot(2,5,ifield-3)
    if iquad==1
        pltytc=loglog(lfine,bly,'r'); hold on
    else
        loglog(lfine,bly,'r');
    end
    xlim([1e3,2e6]);
    ylim([1e-6,2]);
    title(dt.name)
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_l$','interpreter','latex','fontsize',18)
    drawnow
    
    subplot(2,5,ifield+2)
    if iquad==1
        semilogx(lfine,bly,'r');hold on
    else
        loglog(lfine,bly,'r');
    end
    xlim([1e3,2e6]);
    ylim([0,1.2]);
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_l$','interpreter','latex','fontsize',18)
    drawnow
    
end
end
savename=strcat(savedir,'bl');
print(savename,'-dpng');close    

%% plot radial profile
figure
setwinsize(gcf,1500,500)

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'yt/inst2/j0_14/',dt.name,'_',quad,'.fits');
    psfytc=fitsread(file);
    
    profile = radial_prof_pk(psfytc,ones(801),401,401,200);
    r=profile.r*0.7;
    pytc=profile.prof./profile.prof(1);
    eytc=profile.err./profile.prof(1);    
    subplot(2,5,ifield-3)    
    if iquad==1
    loglog(r,pytc,'r.-');hold on
    elseif iquad==2
    errorbar(r,pytc,eytc,'r.-');
    else
    loglog(r,pytc,'r.-');hold on
    end
    xlim([0.7,1e2])
    ylim([1e-4,2])
    title(dt.name)
    xlabel('arcsec')
    ylabel('PSF')
    
    subplot(2,5,ifield+2)
    plot(r,pytc,'r.-');hold on
    if iquad==2
    errorbar(r,pytc,eytc,'r.-');
    end
    
    xlim([0.7,20])
    ylim([0,1.1])
    xlabel('arcsec')
    ylabel('PSF')

    drawnow
end
end
savename=strcat(savedir,'profile');
print(savename,'-dpng');close
