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

%%% load signal PS
load(strcat(weightdir,'avgClsig'),'avgClsig');
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
    %%% Kelsall ZL
    zlname = strcat(zldir,'TM1_40030_zlgrad',dt.name,'.fits');
    zl = fliplr(fitsread_orient(zlname));
    zlmean= mean(zl(:));
    zlgrad_true= zl-zlmean;
    %%% RN 2DCl for RN realization
    Cl2d_ave=darkstat(ifield).Cl2df_ave;
    %%% mask
    maskdir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
    load(strcat(maskdir,'maskin'),'maskin');
    %%% flgiht mean for ph realization
    [flightfr] = get_data_frames...
            (inst,dt.name,'flight',flight,'verbose',0);
    flightfr=flightfr(3:end,:,:);
    slopemap=linfit_map(flightfr,'verbose',0);
    meanmap=mean(slopemap(find(maskin))).*cal;
    clear flightfr slopemap
    %%% mkk
    Mkk=fwfulldat(ifield).Mkk;
    %%% sim sky signal
    [~,~,~,~,lbin]=get_angular_spec(randn(1024),randn(1024),pixscale);
    %Clin=zeros(29);Clin(9:29)=avgClsig;
    Clin=sigCl_extrap_mz14;
    
    sigmap_true = map_from_power_spec(lbin,Clin,1024,1024,pixscale,1);
    [Cls,~,~,ll]=get_Cl(sigmap_true,maskin,Mkk,pixscale,weight);
    %%% PSF conv
    %psfdir='/Users/ytcheng/ciber/doc/20170320_FFsim/psf/';
    %load(strcat(psfdir,'psf_',dt.name),'psfmap');
    load('/Users/ytcheng/ciber/code/sim/psfSWIRE','psfmap');
    
    sigmap=conv2(sigmap_true,psfmap./sum(psfmap(:)));
    Nout=(size(sigmap,1)-size(sigmap_true,1))/2;
    sigmap=sigmap(Nout+1:Nout+1024,Nout+1:Nout+1024);
    [Clsc,~,~,ll]=get_Cl(sigmap,maskin,Mkk,pixscale,weight);
    
    zlgrad=conv2(zlgrad_true,psfmap./sum(psfmap(:)));
    Nout=(size(zlgrad,1)-size(zlgrad_true,1))/2;
    zlgrad=zlgrad(Nout+1:Nout+1024,Nout+1:Nout+1024);
    %%% get bl 
    %load(strcat('/Volumes/HD1TB/CIBER/data/',num2str(flight),...
    %'/dr/TM',num2str(inst),'_',dt.name,'_dr150206.mat'));
    %logbl=interp1(log10(data.psf.l),log10(data.psf.bl),log10(ll),'spline');
    %bl=10.^logbl;
    load('/Users/ytcheng/ciber/code/sim/blSWIRE','blSWIRE');
    bl=blSWIRE(9:29);
    
    siminfo(ifield).name=dt.name;
    siminfo(ifield).nfr=dt.nfr;
    siminfo(ifield).weight=weight;
    siminfo(ifield).zlgrad_true=zlgrad_true;
    siminfo(ifield).zlgrad=zlgrad;
    siminfo(ifield).zlmean=zlmean;
    siminfo(ifield).Cl2d_ave=Cl2d_ave;
    siminfo(ifield).meanmap=meanmap;
    siminfo(ifield).bigmask=maskin;
    siminfo(ifield).Mkk=Mkk;
    siminfo(ifield).sigmap_true=sigmap_true;
    siminfo(ifield).Cls=Cls;
    siminfo(ifield).sigmap=sigmap;
    siminfo(ifield).Clsc=Clsc;    
    siminfo(ifield).bl=bl;
    siminfo(ifield).l=ll;
end
%% 100 runs
tic
nsim=100;
Cl1_arr=zeros(8,nsim,21);
Clx_arr=zeros(size(xmkk,2),nsim,21);
for isim=1:nsim
%%%%% get sim maps %%%%%%
for ifield =[8,7,6,5,4,2,1]
    sigmap=siminfo(ifield).sigmap;
    zlgrad=siminfo(ifield).zlgrad;
    meanmap=siminfo(ifield).meanmap;
    slopemap=(sigmap+zlgrad+meanmap).*labflat./cal;
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
    zlgrad=siminfo(ifield).zlgrad;
    sigmap=siminfo(ifield).sigmap;
    meanmap=siminfo(ifield).meanmap;
    phmap=mapdat(ifield).phmap;
    rnmap=mapdat(ifield).rnmap;

    obs=((sigmap+zlgrad+meanmap).*labflat)+phmap+rnmap;
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
    plane=plane_fit(map1,bigmask);
    map1=(map1-plane).*bigmask;
    
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
    [Cl1,~,~,l]=get_Cl(map1,bigmask,Mkk,pixscale,weight);
    Cl1_arr(ifield,isim,:)=Cl1;
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
simsigCldat.l=l;
meanCl1=squeeze(mean(Cl1_arr(ifield,:,:)));
stdCl1=squeeze(std(Cl1_arr(ifield,:,:)));
simsigCldat.auto(ifield).Cl1_arr=squeeze(Cl1_arr(ifield,:,:));
simCldat.auto(ifield).meanCld=meanCl1';
simCldat.auto(ifield).stdCld=stdCl1';
simsigCldat.siminfo(ifield).name=siminfo(ifield).name;
simsigCldat.siminfo(ifield).Cls=siminfo(ifield).Cls;
simsigCldat.siminfo(ifield).Clsc=siminfo(ifield).Clsc;
simsigCldat.siminfo(ifield).bl=siminfo(ifield).bl;
end

count=0;
for ifield=[8,7,6,5,4,2,1]
    for jfield=[8,7,6,5,4,2,1]
        if jfield<ifield
            count=count+1;
            meanClx=squeeze(mean(Clx_arr(count,:,:)));
            stdClx=squeeze(std(Clx_arr(count,:,:)));
            simsigCldat.cross(count).field=xmkk(count).field;
            simsigCldat.cross(count).Clx_arr=squeeze(Clx_arr(count,:,:));
            simsigCldat.cross(count).meanClx=meanClx';
            simsigCldat.cross(count).stdClx=stdClx';
        end
    end
end
save(strcat(savedir,'simsigCldat'),'simsigCldat');
