inst=1;
pixscale=7;
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1=-3.6;G2=cal./G1;
savedir='/Users/ytcheng/ciber/doc/20170320_FFsim/';
zldir = '/Users/ytcheng/ciber/data/TM1_zl_pmk/';

%%% get Fourier weight
weightdir='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
load(strcat(weightdir,'fwfulldat'),'fwfulldat');

%%% load cross mkk
xmkkdir='/Users/ytcheng/ciber/doc/20170325_alldat/TM1/';
load(strcat(xmkkdir,'xmkk'),'xmkk');
%%% load labflat                
drdir=strcat('/Volumes/HD1TB/CIBER/data/',num2str(flight),'/dr/');
dataname = strcat(drdir,'TM',num2str(inst),'_SWIRE_dr150206.mat');
% lab flat are the same for all fields
load(dataname);labflat=data.flat.labflat;
clear data

%%% load darkstat for RN realization
load(strcat(savedir,'darkstat'),'darkstat');
%% simulation input info

for ifield=[8,7,6,5,4,2,1]
    dt=get_dark_times(flight,inst,ifield);
    %%% Fourier weight
    std_noise=fwfulldat(ifield).std_fCl2d;
    weight=(fftshift(fftshift(1./std_noise)))';
    %%% RN 2DCl for RN realization
    Cl2d_ave=darkstat(ifield).Cl2df_ave;
    %%% mask
    maskdir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
    load(strcat(maskdir,'maskin'),'maskin');

    %%% flgiht mean for ph realization
    load(strcat(maskdir,'flightmap'),'flightmap');
    slopemap=flightmap.rawmapf;
    meanmap=mean(slopemap(find(maskin))).*cal;
    clear flightfr slopemap
    %%% mkk
    Mkk=fwfulldat(ifield).Mkk;

    siminfo(ifield).name=dt.name;
    siminfo(ifield).nfr=dt.nfr;
    siminfo(ifield).weight=weight;
    siminfo(ifield).Cl2d_ave=Cl2d_ave;
    siminfo(ifield).meanmap=meanmap;
    siminfo(ifield).bigmask=maskin;
    siminfo(ifield).Mkk=Mkk;
end
%% 100 runs
tic
nsim=100;
Cl1_arr=zeros(8,nsim,21);
Cln_arr=zeros(8,nsim,21);
Cld_arr=zeros(8,nsim,21);
Clx_arr=zeros(size(xmkk,2),nsim,21);
for isim=1:nsim
    
%%%%% get sim maps %%%%%%
for ifield =[8,7,6,5,4,2,1]
    meanmap=siminfo(ifield).meanmap;
    slopemap=(meanmap).*labflat./cal;
    [phmap]=photonnoise_realization...
        (slopemap,G1,siminfo(ifield).nfr,frate);
    rnmap=readnoise_realization...
                (siminfo(ifield).Cl2d_ave,pixscale,'norand',1);
    mapdat(ifield).phmap=phmap.*cal;
    mapdat(ifield).rnmap=rnmap.*cal;
end

%%%%% get obs map %%%%%    
for ifield=[8,7,6,5,4,2,1]
    bigmask=siminfo(ifield).bigmask;
    meanmap=siminfo(ifield).meanmap;
    phmap=mapdat(ifield).phmap;
    rnmap=mapdat(ifield).rnmap;

    obs=(meanmap.*labflat)+phmap+rnmap;
    obs=obs.*bigmask;
    mapdat(ifield).obsmap=obs;
end
%%%% stack FF %%%%%
for ifield=[8,7,6,5,4,2,1]
FF=zeros(size(labflat));stack_mask=zeros(size(labflat));
FFuse=[0,0,0,1,1,1,1,1];FFuse(ifield)=0;
for jfield=[8,7,6,5,4,2,1]
    if FFuse(jfield)==1
    obs=mapdat(jfield).obsmap;
    bigmask=siminfo(jfield).bigmask;
    mean_obs=mean(obs(find(obs)));
    FF=FF+(obs./sqrt(mean_obs));
    stack_mask=stack_mask+bigmask.*sqrt(mean_obs);
    end
end
FF=FF./stack_mask;FF((find(FF~=FF)))=0;
mapdat(ifield).FF=FF;
end

%%%%% get FF corrected map %%%%%
for ifield=[8,7,6,5,4,2,1]
    bigmask=siminfo(ifield).bigmask;
    phmap=mapdat(ifield).phmap;
    rnmap=mapdat(ifield).rnmap;
    FF=mapdat(ifield).FF;
    
    frnmap=rnmap./FF;
    fphmap=phmap./FF;
    obs=mapdat(ifield).obsmap;
    map1=obs./FF;
    map1(find(map1~=map1))=0;map1(find(map1==inf))=0;
    map1=map1.*bigmask;
    
    mapdat(ifield).frnmap=frnmap;
    mapdat(ifield).fphmap=fphmap;
    mapdat(ifield).map1=map1;
end

%%%%%% get auto spectra %%%%%%
for ifield=[8,7,6,5,4,2,1]
    bigmask=siminfo(ifield).bigmask;
    weight=siminfo(ifield).weight;
    Mkk=siminfo(ifield).Mkk;    
    frnmap=mapdat(ifield).frnmap;
    fphmap=mapdat(ifield).fphmap;
    map1=mapdat(ifield).map1;
    [Cln]=get_Cl(fphmap+frnmap,bigmask,Mkk,pixscale,weight);
    [Cl1,~,~,l]=get_Cl(map1,bigmask,Mkk,pixscale,weight);
    Cld=Cl1-Cln;
    Cl1_arr(ifield,isim,:)=Cl1;
    Cln_arr(ifield,isim,:)=Cln; 
    Cld_arr(ifield,isim,:)=Cld;
end

%%%%%% get cross spectra %%%%%%
count=0;
for ifield=[8,7,6,5,4,2,1]
    bigmaski=siminfo(ifield).bigmask;
    weighti=siminfo(ifield).weight;
    map1i=mapdat(ifield).map1;
    for jfield=[8,7,6,5,4,2,1]
        if jfield<ifield
            count=count+1;
            bigmaskj=siminfo(jfield).bigmask;
            bigmask=bigmaski.*bigmaskj;
            weightj=siminfo(jfield).weight;
            weight=(weighti+weightj)./2;
            Mkk=xmkk(count).Mkk;
            map1j=mapdat(jfield).map1;
            [Clx]=get_Clx(map1i,map1j,bigmask,Mkk,pixscale,weight);
            Clx_arr(count,isim,:)=Clx;
        end
        
    end
    
end
disp(sprintf('isim=%d,t=%.1f min',isim,toc/60));
end
%%
for ifield=[8,7,6,5,4,2,1]
simCldat.l=l;
meanCld=squeeze(mean(Cld_arr(ifield,:,:)));
stdCld=squeeze(std(Cld_arr(ifield,:,:)));
simCldat.auto(ifield).meanCld=meanCld';
simCldat.auto(ifield).stdCld=stdCld';
end

count=0;
for ifield=[8,7,6,5,4,2,1]
    for jfield=[8,7,6,5,4,2,1]
        if jfield<ifield
            count=count+1;
            meanClx=squeeze(mean(Clx_arr(count,:,:)));
            stdClx=squeeze(std(Clx_arr(count,:,:)));
            simCldat.cross(count).field=xmkk(count).field';
            simCldat.cross(count).meanClx=meanClx';
            simCldat.cross(count).stdClx=stdClx';
        end
    end
end
save(strcat(savedir,'simCldat'),'simCldat');
%% plot
for ifield =[8,7,6,5,4,2,1]
Cl1=squeeze(Cl1_arr(ifield,:,:));
Cln=squeeze(Cln_arr(ifield,:,:));
l=simCldat.l;

fig=figure;                
pltCl1=errorbar(l,l.*(l+1).*prctile(Cl1,50)./2./pi,...
     l.*(l+1).*(prctile(Cl1,50)-prctile(Cl1,16))./2./pi,...
     l.*(l+1).*(prctile(Cl1,84)-prctile(Cl1,50))./2./pi,...
                    'or','markersize',3,'DisplayName',...
'$M\frac{FF}{\hat{FF}}+\frac{\delta RN}{\hat{FF}}+\frac{\delta ph}{\hat{FF}}$');
hold on

y1=(l.*(l+1).*(prctile(Cln,16))./2./pi);
y2=(l.*(l+1).*(prctile(Cln,84))./2./pi);
pltCln=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none','DisplayName',...
'$\frac{\delta RN}{\hat{FF}}+\frac{\delta ph}{\hat{FF}}$');hold on


xlim([1e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','southeast');
set(h,'Interpreter','latex','FontSize',20);
legend boxoff
title(siminfo(ifield).name); 
end
