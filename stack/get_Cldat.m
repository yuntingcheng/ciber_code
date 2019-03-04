function get_Cldat(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the flight full integration PS, noise PS. Used for fit RN and ph
% level with and bl diagnostic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);
pixscale=7;
cp=get_cal_params('flight',flight);
frate=cp(inst).framerate;

loaddir=sprintf('%sTM%d/',mypaths.filtmap,inst);
load(sprintf('%s/darklongdat',loaddir),'darklongdat');
DCtemplate=darklongdat.DCtemplate; clear darklongdat

loaddir1=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sfwdat',loaddir1),'fwdat');
load(sprintf('%sFFdat',loaddir1),'FFdat');
load(sprintf('%smaskdat',loaddir1),'maskdat');
load(sprintf('%smkkdat',loaddir1),'mkkdat');
load(sprintf('%sbldat',loaddir1),'bldat');
load(sprintf('%ssimCldat',loaddir1),'simCldat');

savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/blfit/TM',...
    num2str(inst),'/');
%% get the flight, RN spectrum
for ifield=4:8
    disp(sprintf('ifield=%d',ifield));
    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'flightmap'),'flightmap');
    dt=get_dark_times(flight,inst,ifield);

    fw=fwdat(ifield).fw_filt;
    mkk=mkkdat.auto(ifield).mkk_wfull;
    bigmask=maskdat.mask(ifield).bigmask;
    FF=FFdat(ifield).FF;

    %%% flight PS
    flightf=flightmap.filtmapf;
    flightf=(flightf-DCtemplate)./FF.*bigmask;
    plane=plane_fit(flightf,bigmask);
    flightf=flightf-plane.*bigmask;
    flightf(find(flightf~=flightf))=0;
    flightf(find(flightf==-inf))=0;
    flightf(find(flightf==inf))=0;
    calmap=flightf;
    [flightCl,~,~,l]=get_Cl(calmap,bigmask,mkk,pixscale,fw);
    
    Cldat(ifield).l=l;
    Cldat(ifield).flightCl=flightCl;
        
    %%% FF bias
    Cldat(ifield).rnClFF=simCldat.nbias.auto(ifield).rnavgCl;
    Cldat(ifield).drnClFF=simCldat.nbias.auto(ifield).rnstdCl;
    Cldat(ifield).phClFF=simCldat.nbias.auto(ifield).phavgCl;
    Cldat(ifield).dphClFF=simCldat.nbias.auto(ifield).phstdCl;

    %%% read noise PS
    rnCl_arr=zeros(numel(dt.time),21);
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');
        filtmapf=labmap.filtmapf;
        filtmapf=(filtmapf-DCtemplate)./FF;
        [~,mask]=get_skymap(filtmapf,bigmask,4,5);
        [Cl]=get_Cl(filtmapf,mask,mkk,pixscale,fw);
        rnCl_arr(i,:)=Cl;    
    end
    Cldat(ifield).rnCl=squeeze(mean(rnCl_arr));
    Cldat(ifield).drnCl=squeeze(std(rnCl_arr));

    %%% photon noise PS with G1=-1
    slopemap=flightmap.rawmapf;
    slopemean=mean(slopemap(find(bigmask)));
    phCl_arr=zeros(100,21);
    for i=1:100
        phmap=photonnoise_realization(slopemean.*ones(1024),-1,dt.nfr,frate);
        [phCl]=get_Cl(phmap./FF,bigmask,mkk,pixscale,fw);
        phCl_arr(i,:)=phCl;
    end
    Cldat(ifield).phCl=squeeze(mean(phCl_arr));
    Cldat(ifield).dphCl=squeeze(std(phCl_arr));
    
    
    %%% get bl mz
    load(strcat(mypaths.release,'TM',num2str(inst),'_',dt.name,'_dr150206.mat'));
    lfine=data.psf.l;
    blfine=data.psf.bl;
    lfine=lfine(find(blfine==blfine));
    blfine=blfine(find(blfine==blfine));
    logbl=interp1(log10(lfine),log10(blfine),log10(l),'pchip','extrap');
    blmz=10.^logbl;
    Cldat(ifield).blmz=blmz;
    
end
save(strcat(savedir,'Cldat'),'Cldat');
return