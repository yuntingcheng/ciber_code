%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (TM1) Take the PSF map from 2MASS stacking,
% plot the PSF map, Bl, angular average radial profile.
% Comparing with the PSF from PMK. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
loaddir='/Users/ytcheng/ciber/doc/20170617_Stacking/psf/';
loaddir1='/Users/ytcheng/ciber/data/iband_psf_quads/';
savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/TM1/';
quad_arr=['A','B','C','D'];
%% plot the psf
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
figure
setwinsize(gcf,1200,750)
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'yt/inst1/j0_14/',dt.name,'_',quad,'.fits');
    psfytc=fitsread(file);
    file=strcat(loaddir,'pmk/',dt.name,'_',quad,'.fits');
    psfpmk=fitsread(file);
    psfpmk1=fitsread(strcat(loaddir1,dt.name,'_',quad,'.fits'));
    
    subplot(3,4,iquad)
    imageclip(psfpmk1(300:500,300:500));
    v=caxis;
    title(strcat(dt.name,'\_',quad))
    subplot(3,4,iquad+4)
    imageclip(psfpmk(300:500,300:500));
    caxis(v)
    subplot(3,4,iquad+8)
    imageclip(psfytc(300:500,300:500));
    %caxis(v)
    
end
savename=strcat(savedir,'psfmap_',dt.name);
print(savename,'-dpng');close
end
%% plot bl
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
figure
setwinsize(gcf,1000,400)
for iquad=1:4
    
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'yt/inst1/j0_14/',dt.name,'_',quad,'.fits');
    psfytc=fitsread(file);
    file=strcat(loaddir,'pmk/',dt.name,'_',quad,'.fits');
    psfpmk=fitsread(file);
    psfpmk1=fitsread(strcat(loaddir1,dt.name,'_',quad,'.fits'));
    
    psfytc(find(psfytc~=psfytc))=0;
    psfpmk(find(psfpmk~=psfpmk))=0;
    psfpmk1(find(psfpmk1~=psfpmk1))=0;
    
    [bly,lfine]=get_angular_spec(psfytc,psfytc,0.7,'nbins',30);
    [blp]=get_angular_spec(psfpmk,psfpmk,0.7,'nbins',30);
    [blp1]=get_angular_spec(psfpmk1,psfpmk1,0.7,'nbins',30);

    lfine=lfine(find(bly));
    bly=bly(find(bly));blp=blp(find(blp));blp1=blp1(find(blp1));
    bly=bly./mean(bly(1:3));
    blp=blp./mean(blp(1:3));
    blp1=blp1./mean(blp1(1:3));
    
    subplot(1,2,1)
    if iquad==1
        pltytc=loglog(lfine,bly,'r');
        hold on
        pltpmk=loglog(lfine,blp,'g');
        pltpmk1=loglog(lfine,blp1,'b');
    else
        loglog(lfine,bly,'r');
        hold on
        loglog(lfine,blp,'g');
        loglog(lfine,blp1,'b');
    end
    xlim([1e3,2e6]);
    ylim([1e-6,2]);
    title(dt.name)
    legend([pltytc,pltpmk,pltpmk1],...
    {'new psf','old psf code','old psf'},'location','southwest');
    legend boxoff
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_l$','interpreter','latex','fontsize',18)
    drawnow
    
    subplot(1,2,2)
    if iquad==1
        semilogx(lfine,bly,'r');
        hold on
        loglog(lfine,blp,'g');
        loglog(lfine,blp1,'b');
    else
        loglog(lfine,bly,'r');
        hold on
        loglog(lfine,blp,'g');
        loglog(lfine,blp1,'b');
    end
    xlim([1e3,2e6]);
    ylim([0,1.2]);
    title(dt.name);
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_l$','interpreter','latex','fontsize',18)
    drawnow
    
end
savename=strcat(savedir,'bl_',dt.name);
print(savename,'-dpng');%close    

end

%% plot radial profile
figure
setwinsize(gcf,1500,500)

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'yt/inst1/j0_14/',dt.name,'_',quad,'.fits');
    psfytc=fitsread(file);
    psfpmk1=fitsread(strcat(loaddir1,dt.name,'_',quad,'.fits'));
    
    profile = radial_prof_pk(psfytc,ones(801),401,401,200);
    r=profile.r*0.7;
    pytc=profile.prof./profile.prof(1);
    eytc=profile.err./profile.prof(1);
    
    profile = radial_prof_pk(psfpmk1,ones(801),401,401,200);
    ppmk=profile.prof./profile.prof(1);
    epmk=profile.err./profile.prof(1);
    
    subplot(2,5,ifield-3)    
    if iquad==1
    loglog(r,pytc,'r.-','DisplayName','PSF new');hold on
    loglog(r,ppmk,'b.-','DisplayName','PSF old');
    elseif iquad==2
    errorbar(r,pytc,eytc,'r.-');
    errorbar(r,ppmk,epmk,'b.-');
    else
    loglog(r,pytc,'r.-');hold on
    loglog(r,ppmk,'b.-');        
    end
    xlim([0.7,1e2])
    ylim([1e-4,2])
    title(dt.name)
    xlabel('arcsec')
    ylabel('PSF')
    h=legend('show','Location','southwest');
    set(h,'fontsize',10)
    legend boxoff

    
    subplot(2,5,ifield+2)
    plot(r,pytc,'r.-');hold on
    plot(r,ppmk,'b.-');
    if iquad==2
    errorbar(r,pytc,eytc,'r.-');
    errorbar(r,ppmk,epmk,'b.-');    
    end
    
    xlim([0.7,20])
    ylim([0,1.1])
    xlabel('arcsec')
    ylabel('PSF')

end
end

savename=strcat(savedir,'profile');
print(savename,'-dpng');%close
