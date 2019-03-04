%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ihl stacking maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
ifield=8;
npix=800;


savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/UKIDSS/'));
quad_arr=['A','B','C','D'];
%% get the PSF
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');

load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

bestparam=fitpsfdat(ifield).bestparam_norm;
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
%% plot the stacking map
for m=14:20

dt=get_dark_times(flight,inst,ifield);

for itype=1:2
    if itype==1
        type='ciber';
    else
        type='ukidss';
    end
      loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/',char(type),'/'));
  
stackstot=zeros(npix*2+1);
stackgtot=zeros(npix*2+1);
mstackstot=zeros(npix*2+1);
mstackgtot=zeros(npix*2+1);
for iquad=1:4
    quad=quad_arr(iquad);
    
    stacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
        '_stampers',num2str(m),'.fits'));
    stackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
        '_stamperg',num2str(m),'.fits'));
    mstacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
        '_maskstampers',num2str(m),'.fits'));
    mstackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
        '_maskstamperg',num2str(m),'.fits'));
    
    stackstot=stackstot+stacks;
    stackgtot=stackgtot+stackg;
    mstackstot=mstackstot+mstacks;
    mstackgtot=mstackgtot+mstackg;
    
    
    
end
    if itype==1
        stackscb=stackstot./mstackstot;
        stackgcb=stackgtot./mstackgtot;
    else
        stacksuk=stackstot./mstackstot;
        stackguk=stackgtot./mstackgtot;
    end

end
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(stackscb(npix-100:npix+100,npix-100:npix+100));
title('CIBER stars');
subplot(2,2,2)
imageclip(stackgcb(npix-100:npix+100,npix-100:npix+100));
title('CIBER gals');
subplot(2,2,3)
imageclip(stacksuk(npix-100:npix+100,npix-100:npix+100));
title('UKIDSS stars');
subplot(2,2,4)
imageclip(stackguk(npix-100:npix+100,npix-100:npix+100));
title('UKIDSS gals');

savename=strcat(savedir,dt.name,'_stackmaps_allquad',num2str(m));
print(savename,'-dpng');%close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot radial profile
figure
setwinsize(gcf,1000,800)
for itype=1:2
    if itype==1
        stacks=stackscb;
        stackg=stackgcb;
    else
        stacks=stacksuk;
        stackg=stackguk;
    end
    

    norm = stacks(npix+1,npix+1);  
    profile = radial_prof(stacks,ones(2*npix+1),npix+1,npix+1,1,25);
    r_arr=profile.r*0.7;
    profs_arr=(profile.prof)./norm;
    errs_arr=profile.err./norm;

    norm = stackg(npix+1,npix+1);  
    profile = radial_prof(stackg,ones(2*npix+1),npix+1,npix+1,1,25);
    profg_arr=(profile.prof)./norm;
    errg_arr=profile.err./norm;


    if itype==1
        profscb_arr=profs_arr;
        errscb_arr=errs_arr;
        profgcb_arr=profg_arr;
        errgcb_arr=errg_arr;
    else
        profsuk_arr=profs_arr;
        errsuk_arr=errs_arr;
        profguk_arr=profg_arr;
        errguk_arr=errg_arr;
    end

    subplot(2,2,1)
    if itype==1
        errorbar(r_arr,profs_arr,errs_arr,'r.-','DisplayName','CBstars');hold on
        errorbar(r_arr,profg_arr,errg_arr,'b.-','DisplayName','CBgals');
        else
        errorbar(r_arr,profs_arr,errs_arr,'m.-','DisplayName','UKstars');
        errorbar(r_arr,profg_arr,errg_arr,'c.-','DisplayName','UKgals');     
    end
    
    
    subplot(2,2,2)
    if itype==1
        loglog(r_arr,profs_arr,'r.-');hold on
        loglog(r_arr,profg_arr,'b.-');
        errorbar(r_arr,profs_arr,errs_arr,'r.-');
        errorbar(r_arr,profg_arr,errg_arr,'b.-');
    else
        loglog(r_arr,profs_arr,'m.-');hold on
        loglog(r_arr,profg_arr,'c.-');
        errorbar(r_arr,profs_arr,errs_arr,'m.-');
        errorbar(r_arr,profg_arr,errg_arr,'c.-');        
    end

end

subplot(2,2,1)
title(strcat(dt.name,',',num2str(m),'<m<',num2str(m+1)))
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

subplot(2,2,2)
loglog(r_arr,profpsf_arr,'k--');
xlim([4e-1,7e2])
ylim([1e-4,2e0])
xlabel('arcsec')
ylabel('<I_{stack}>')


subplot(2,2,3)
errorbar(r_arr.*0.98,profsuk_arr,errsuk_arr,'r.-','DisplayName','2Mstars');hold on
errorbar(r_arr.*1.02,profguk_arr,errguk_arr,'b.-','DisplayName','2Mgals');
loglog(r_arr,profpsf_arr,'k--','DisplayName','PSF');
set(gca, 'XScale', 'log','YScale', 'log');
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
xlim([5e-1,1e2])
ylim([5e-2,2e0])

title('UKIDSS')
xlabel('arcsec')
ylabel('<I_{stack}>')

subplot(2,2,4)
errorbar(r_arr.*0.98,profscb_arr,errscb_arr,'r.-','DisplayName','CBstars');hold on
errorbar(r_arr.*1.02,profgcb_arr,errgcb_arr,'b.-','DisplayName','CBgals');
loglog(r_arr,profpsf_arr,'k--','DisplayName','PSF');
set(gca, 'XScale', 'log','YScale', 'log');
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
xlim([5e-1,1e2])
ylim([5e-2,2e0])
title('CIBER')
xlabel('arcsec')
ylabel('<I_{stack}>')

savename=strcat(savedir,dt.name,'_rprof_allquad',num2str(m));
print(savename,'-dpng');%close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot radial profile offset
figure
setwinsize(gcf,1000,800)
for itype=1:2
    if itype==1
        stacks=stackscb;
        stackg=stackgcb;
    else
        stacks=stacksuk;
        stackg=stackguk;
    end
    

    norm = stacks(npix+1,npix+1);  
    profile = radial_prof(stacks,ones(2*npix+1),npix+1,npix+1,1,25);
    r_arr=profile.r*0.7;
    sqrtn_arr=sqrt(profile.npix);
    off=mean(profile.prof(end-5:end));
    profs_arr=profile.prof-off;
    norm=norm-off;
    profs_arr=(profs_arr)./norm;
    errs_arr=profile.err./norm;

    norm = stackg(npix+1,npix+1);
    profile = radial_prof(stackg,ones(2*npix+1),npix+1,npix+1,1,25);
    sqrtn_arr=sqrt(profile.npix);
    off=mean(profile.prof(end-5:end));
    profg_arr=profile.prof-off;
    norm=norm-off;
    profg_arr=(profg_arr)./norm;
    errg_arr=profile.err./norm;


    if itype==1
        profscb_arr=profs_arr;
        errscb_arr=errs_arr;
        profgcb_arr=profg_arr;
        errgcb_arr=errg_arr;
    else
        profsuk_arr=profs_arr;
        errsuk_arr=errs_arr;
        profguk_arr=profg_arr;
        errguk_arr=errg_arr;
    end

    subplot(2,2,1)
    if itype==1
        errorbar(r_arr,profs_arr,errs_arr,'r.-','DisplayName','CBstars');hold on
        errorbar(r_arr,profg_arr,errg_arr,'b.-','DisplayName','CBgals');
        else
        errorbar(r_arr,profs_arr,errs_arr,'m.-','DisplayName','UKstars');
        errorbar(r_arr,profg_arr,errg_arr,'c.-','DisplayName','UKgals');     
    end
    
    
    subplot(2,2,2)
    if itype==1
        loglog(r_arr,profs_arr,'r.-');hold on
        loglog(r_arr,profg_arr,'b.-');
        errorbar(r_arr,profs_arr,errs_arr,'r.-');
        errorbar(r_arr,profg_arr,errg_arr,'b.-');
    else
        loglog(r_arr,profs_arr,'m.-');hold on
        loglog(r_arr,profg_arr,'c.-');
        errorbar(r_arr,profs_arr,errs_arr,'m.-');
        errorbar(r_arr,profg_arr,errg_arr,'c.-');        
    end

end

subplot(2,2,1)
title(strcat(dt.name,',',num2str(m),'<m<',num2str(m+1)))
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

subplot(2,2,2)
loglog(r_arr,profpsf_arr,'k--');
xlim([4e-1,7e2])
ylim([1e-4,2e0])
xlabel('arcsec')
ylabel('<I_{stack}>')


subplot(2,2,3)
errorbar(r_arr.*0.98,profsuk_arr,errsuk_arr,'r.-','DisplayName','CBstars');hold on
errorbar(r_arr.*1.02,profguk_arr,errguk_arr,'b.-','DisplayName','CBgals');
loglog(r_arr,profpsf_arr,'k--','DisplayName','PSF');
set(gca, 'XScale', 'log','YScale', 'log');
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
xlim([5e-1,1e2])
ylim([5e-2,2e0])

title('UKIDSS')
xlabel('arcsec')
ylabel('<I_{stack}>')

subplot(2,2,4)
errorbar(r_arr.*0.98,profscb_arr,errscb_arr,'r.-','DisplayName','UKstars');hold on
errorbar(r_arr.*1.02,profgcb_arr,errgcb_arr,'b.-','DisplayName','UKgals');
loglog(r_arr,profpsf_arr,'k--','DisplayName','PSF');
set(gca, 'XScale', 'log','YScale', 'log');
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
xlim([5e-1,1e2])
ylim([5e-2,2e0])
title('CIBER')
xlabel('arcsec')
ylabel('<I_{stack}>')


savename=strcat(savedir,dt.name,'_rprof_allquadoff',num2str(m));
print(savename,'-dpng');%close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% all stars & gals
figure
setwinsize(gcf,1000,400)

for isrc=[1,2]
    if isrc==1
        name='stars';
    else
        name='galaxies';
    end
    
subplot(1,2,isrc);
for icase=[1,2]
if icase==1
loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/ukidss/')); 
c='b';
else
loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/CIBER/')); 
c='r';
end

for m=16:20
stackstot=zeros(npix*2+1);
mstackstot=zeros(npix*2+1);
for iquad=1:4
    quad=quad_arr(iquad);
    
    if isrc==1
        stacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_stampers',num2str(m),'.fits'));
        mstacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_maskstampers',num2str(m),'.fits'));
    else
        stacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_stamperg',num2str(m),'.fits'));
        mstacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_maskstamperg',num2str(m),'.fits'));
    end
    
        
    stackstot=stackstot+stacks;
    mstackstot=mstackstot+mstacks;
      
end
stackscb=stackstot./mstackstot;
stacks=stackscb;
    

norm = stacks(npix+1,npix+1);
profile = radial_prof(stacks,ones(2*npix+1),npix+1,npix+1,1,25);
r_arr=profile.r*0.7;
off=mean(profile.prof(end-5:end));
profs_arr=profile.prof-off;
norm=norm-off;
profs_arr=(profs_arr)./norm;
errs_arr=profile.err./norm;


profscb_arr=profs_arr;
errscb_arr=errs_arr;

norms=(profscb_arr(1));
loglog(r_arr,(profscb_arr)./norms,'-','color',c);hold on
errorbar(r_arr,(profscb_arr)./norms,errscb_arr./norms,'.','color',c);
drawnow
end

end
loglog(r_arr,profpsf_arr,'k--','linewidth',3);
xlim([5e-1,1e3])
ylim([1e-4,2e0])
title(name)
xlabel('arcsec')
ylabel('<I_{stack}>')
end
savename=strcat(savedir,dt.name,'_rprof_srcs');
print(savename,'-dpng');%close

