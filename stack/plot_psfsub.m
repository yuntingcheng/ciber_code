%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (TM1/TM2) Consistency check:Comparing PSF map, radial profile, and Bl
% from whole samples (0<Jmag<14) with bright and faint subsamples
% (split by Jmag=13)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flight=40030;
inst=1;
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf/yt/inst',...
    num2str(inst),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/plots/TM',...
        num2str(inst),'/subsample/');
quad_arr=['A','B','C','D'];
%% plot the psf
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
figure
setwinsize(gcf,1200,750)
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'j0_14/',dt.name,'_',quad,'.fits');
    psf=fitsread(file);
    file=strcat(loaddir,'j0_13/',dt.name,'_',quad,'.fits');
    psf1=fitsread(file);
    file=strcat(loaddir,'j13_14/',dt.name,'_',quad,'.fits');
    psf2=fitsread(file);
    
    subplot(3,4,iquad)
    imageclip(psf(300:500,300:500));
    v=caxis;
    title(strcat(dt.name,'\_',quad))
    subplot(3,4,iquad+4)
    imageclip(psf1(300:500,300:500));
    caxis(v)
    title('0<j<13')
    subplot(3,4,iquad+8)
    imageclip(psf2(300:500,300:500));
    caxis(v)
    title('13<j<14')
    
end
savename=strcat(savedir,'psfmap_',dt.name);
print(savename,'-dpng');close
end
%% plot radial profile
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
figure
setwinsize(gcf,1200,500)
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'j0_14/',dt.name,'_',quad,'.fits');
    psf=fitsread(file);
    file=strcat(loaddir,'j0_13/',dt.name,'_',quad,'.fits');
    psf1=fitsread(file);
    file=strcat(loaddir,'j13_14/',dt.name,'_',quad,'.fits');
    psf2=fitsread(file);
    
    profile = radial_prof_pk(psf,ones(801),401,401,200);
    r=profile.r*0.7;
    p=profile.prof./profile.prof(1);
    e=profile.err./profile.prof(1);
    
    profile = radial_prof_pk(psf1,ones(801),401,401,200);
    p1=profile.prof./profile.prof(1);
    e1=profile.err./profile.prof(1);

    profile = radial_prof_pk(psf2,ones(801),401,401,200);
    p2=profile.prof./profile.prof(1);
    e2=profile.err./profile.prof(1);

    psf(find(psf~=psf))=0;
    psf1(find(psf1~=psf1))=0;
    psf2(find(psf2~=psf2))=0;
    
    subplot(2,4,iquad)
    loglog(r,p,'k.-');hold on
    loglog(r,p1,'b.-');hold on
    loglog(r,p2,'r.-');hold on
    xlim([0.7,1e2])
    ylim([1e-4,2])
    title(strcat(dt.name,'\_',quad))
    xlabel('arcsec')
    ylabel('PSF')
    legend({'0<J<14','0<J<13','13<J<14'},'location','northeast');
    legend boxoff
    
    subplot(2,4,iquad+4)
    plot(r,p,'k.-');hold on
    plot(r,p1,'b.-');
    plot(r,p2,'r.-');
    xlim([0.7,20])
    ylim([0,1.1])
    xlabel('arcsec')
    ylabel('PSF')
    legend({'0<J<14','0<J<13','13<J<14'},'location','northeast');
    legend boxoff
    
    drawnow
end
savename=strcat(savedir,'profile_',dt.name);
print(savename,'-dpng');close
end
%% plot bl
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
figure
setwinsize(gcf,1200,500)
for iquad=1:4
    quad=quad_arr(iquad);
    
    file=strcat(loaddir,'j0_14/',dt.name,'_',quad,'.fits');
    psf=fitsread(file);
    file=strcat(loaddir,'j0_13/',dt.name,'_',quad,'.fits');
    psf1=fitsread(file);
    file=strcat(loaddir,'j13_14/',dt.name,'_',quad,'.fits');
    psf2=fitsread(file);
    
    psf(find(psf~=psf))=0;
    psf1(find(psf1~=psf1))=0;
    psf2(find(psf2~=psf2))=0;
    
    [blp,lfine]=get_angular_spec(psf,psf,0.7,'nbins',30);
    [blp1]=get_angular_spec(psf1,psf1,0.7,'nbins',30);
    [blp2]=get_angular_spec(psf2,psf2,0.7,'nbins',30);

    lfine=lfine(find(blp));
    blp=blp(find(blp));blp1=blp1(find(blp1));blp2=blp2(find(blp2));
    blp=blp./mean(blp(1:3));
    blp1=blp1./mean(blp1(1:3));
    blp2=blp2./mean(blp2(1:3));

    subplot(2,4,iquad)
    loglog(lfine,blp,'k+');
    hold on
    loglog(lfine,blp1,'bo');
    loglog(lfine,blp2,'ro');
    xlim([1e3,2e6]);
    ylim([1e-6,2]);
    title(strcat(dt.name,'\_',quad))
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_l$','interpreter','latex','fontsize',18)
    legend({'0<J<14','0<J<13','13<J<14'},'location','southwest');
    legend boxoff
    
    subplot(2,4,iquad+4)
    semilogx(lfine,blp,'k+');
    hold on
    plot(lfine,blp1,'bo');
    plot(lfine,blp2,'ro');
    xlim([1e3,2e6]);
    ylim([0,1.2]);
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_l$','interpreter','latex','fontsize',18)
    legend({'0<J<14','0<J<13','13<J<14'},'location','northeast');
    legend boxoff
    drawnow
end
savename=strcat(savedir,'bl_',dt.name);
print(savename,'-dpng');close
end
