%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ihl stacking maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
ifield=8;
npix=800;

dt=get_dark_times(flight,inst,ifield);

savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/2MASS/'));
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
%% get the number counts data
ncountfile=strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/tmass_2m/',dt.name,'_stackcounts.txt');
T=readtable(ncountfile);
mbot_arr=T{:,3}';
counts_arr=T{:,5}';
countg_arr=T{:,6}';
%% plot the stacking map
for m=11:17
counts=sum(counts_arr(find(mbot_arr==m)));
countg=sum(countg_arr(find(mbot_arr==m)));
for itype=1:2
    if itype==1
        type='ciber_2m';
    else
        type='tmass_2m';
    end
      loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/',char(type),'/'));
  
    stackstot=zeros(npix*2+1);
    stackgtot=zeros(npix*2+1);
    mstackstot=zeros(npix*2+1);
    mstackgtot=zeros(npix*2+1);
    for iquad=1:4
        quad=quad_arr(iquad);
        if itype==1
            stacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
                '_stampers',num2str(m),'.fits'));
            stackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
                '_stamperg',num2str(m),'.fits'));
        else
            stacks=fitsread(strcat(loaddir,dt.name,'_',quad,...
                '_stampers',num2str(m),'.fits'));
            stackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
                '_stamperg',num2str(m),'.fits'));            
        end
        
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
        stacks2m=stackstot./mstackstot;
        stackg2m=stackgtot./mstackgtot;
    end
    
end
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(stackscb(npix-100:npix+100,npix-100:npix+100));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,2)
imageclip(stackgcb(npix-100:npix+100,npix-100:npix+100));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
subplot(2,2,3)
imageclip(stacks2m(npix-100:npix+100,npix-100:npix+100));
title(strcat('2MASS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,4)
imageclip(stackg2m(npix-100:npix+100,npix-100:npix+100));
title(strcat('2MASS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

savename=strcat(savedir,dt.name,'_stackmaps_allquad',num2str(m));
print(savename,'-dpng');%close

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot radial profile
figure
setwinsize(gcf,1000,400)
ymin=[];
for itype=1:2
    if itype==1
        stacks=stackscb;
        stackg=stackgcb;
    else
        stacks=stacks2m;
        stackg=stackg2m;
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
    
    %%% force last bin to 0
    %profs_arr=profs_arr-profs_arr(end);
    %profg_arr=profg_arr-profg_arr(end);
    
    ymin=[ymin profs_arr profg_arr];
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
    
    
    
    subplot(1,2,1)
    if itype==1
        errorbar(r_arr,profs_arr,errs_arr,'r.-','DisplayName','CBstars');hold on
        errorbar(r_arr,profg_arr,errg_arr,'b.-','DisplayName','CBgals');
        else
        errorbar(r_arr,profs_arr,errs_arr,'m.-','DisplayName','2Mstars');
        errorbar(r_arr,profg_arr,errg_arr,'c.-','DisplayName','2Mgals');     
    end
    
    
    subplot(1,2,2)
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

subplot(1,2,1)
title(strcat(dt.name,',',num2str(m),'<m<',num2str(m+1)))
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
ylim([-0.2,1.4])
set(gca, 'XScale', 'log')
h=legend('show','Location','southwest');
text(100,1.3,strcat('stack',{' '}, num2str(counts),{' '},'stars'));
text(100,1.2,strcat('stack',{' '}, num2str(countg),{' '},'gals'));

set(h,'fontsize',10)
legend boxoff

ymin=ymin(find(ymin>0));
ymin=floor(log10(min(ymin)));
ymin=10^ymin;
subplot(1,2,2)
loglog(r_arr,profpsf_arr,'k--');
xlim([4e-1,7e2])
ylim([ymin,2e0])
xlabel('arcsec')
ylabel('<I_{stack}>')



savename=strcat(savedir,dt.name,'_rprof_allquad',num2str(m));
%savename=strcat(savedir,dt.name,'_rprof_allquad',num2str(m),'off');

print(savename,'-dpng');%close
end
