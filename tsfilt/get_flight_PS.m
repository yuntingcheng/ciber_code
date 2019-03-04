function get_flight_PS(flight,inst,G1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the flight full integration power spectrum,
% noise power spectrum, and beam bl
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;

loaddir=sprintf('%sTM%d/',mypaths.filtmap,inst);
load(sprintf('%s/darklongdat',loaddir),'darklongdat');
DCtemplate=darklongdat.DCtemplate; clear darklongdat

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/fwdat',savedir),'fwdat');
load(sprintf('%s/FFdat',savedir),'FFdat');
load(sprintf('%s/maskdat',savedir),'maskdat');
load(sprintf('%s/mkkdat',savedir),'mkkdat');
load(sprintf('%s/bldat',savedir),'bldat');
load(sprintf('%s/simCldat',savedir),'simCldat');
%%
for ifield=4:8
    disp(sprintf('ifield=%d',ifield));
    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'flightmap'),'flightmap');
    dt=get_dark_times(flight,inst,ifield);

    fw=fwdat(ifield).fw_filt;
    mkk=mkkdat.auto(ifield).mkk_wfull;
    bigmask=maskdat.mask(ifield).bigmask;

    slopemap=flightmap.rawmapf;
    slopemean=mean(slopemap(find(bigmask)));

    FF=FFdat(ifield).FF;
%     FFcorfac=simCldat.auto(ifield).meanCld;
%     FFcorfacerr=simCldat.auto(ifield).stdCld;
%     FFcorfac=simCldat.nbias.auto(ifield).rnavgCl+...
%                simCldat.nbias.auto(ifield).phavgCl;
%     FFcorfacerr=sqrt(simCldat.nbias.auto(ifield).rnstdCl.^2 + ...
%         simCldat.nbias.auto(ifield).rnstdCl.^2);
    FFcorfac=simCldat.nbias.auto(ifield).navgCl;
    FFcorfacerr=simCldat.nbias.auto(ifield).nstdCl;

    FFcorfac = FFcorfac.*cal.*cal;
    FFcorfacerr = FFcorfacerr.*cal.*cal;
    %%% Noise PS
    nCl_arr=zeros(numel(dt.time),21);
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');
        filtmapf=labmap.filtmapf;
        phmap=photonnoise_realization...
            (slopemean.*ones(1024),G1,dt.nfr,frate);
        filtmapf=(filtmapf-DCtemplate)+phmap;
        [~,mask]=get_skymap(filtmapf,bigmask,4,5);
        plane=plane_fit(filtmapf,mask);
        filtmapf = (filtmapf - plane).*mask;
        [Cl]=get_Cl(filtmapf,mask,mkk,pixscale,fw);
        nCl_arr(i,:)=Cl.*cal.^2;    
    end
    
    %%% flight PS
    flightf=flightmap.filtmapf;
    flightf=((flightf-DCtemplate)./FF).*bigmask;
    plane=plane_fit(flightf,bigmask);
    flightf=flightf-plane.*bigmask;
    flightf(flightf~=flightf)=0;
    flightf(flightf==-inf)=0;
    flightf(flightf==inf)=0;
    calmap=flightf.*cal;
    flightmap.calmap=calmap;
    [Clflight,~,~,l]=get_Cl(calmap,bigmask,mkk,pixscale,fw);
    save(strcat(loaddir,'flightmap'),'flightmap');
    %%% get bl 
%     blfine=bldat.bl(ifield).bl;
%     lfine=bldat.bl(ifield).l;
%     lfine=lfine(find(blfine==blfine));
%     blfine=blfine(find(blfine==blfine));
%     % add a lowest l point to make sure 
%     % interp converge to 1 at low l
%     blfine=[1 blfine];lfine=[l(1) lfine];
%     
%     % no bl data from PMK, use MZ bl
%     if ifield<4
%         load(strcat('/Volumes/HD1TB/CIBER/data/',num2str(flight),...
%         '/dr/TM',num2str(inst),'_',dt.name,'_dr150206.mat'));
%         lfine=data.psf.l;
%         blfine=data.psf.bl;
%         lfine=lfine(find(blfine==blfine));
%         blfine=blfine(find(blfine==blfine));
%     end    

%     logbl=interp1(log10(lfine),log10(blfine),log10(l),'pchip','extrap');
%     bl=10.^logbl;

    bldata = bldat.bl(ifield).bl;
    ldata = bldat.bl(ifield).l;
    bldata(bldata~=bldata) = 1;
    logbl=interp1(log10(ldata),log10(bldata),log10(l),'pchip','extrap');
    bl=10.^logbl;
    
    Cl_blcor=(Clflight-FFcorfac-mean(nCl_arr))./bl;
    Cl_blcor_err=sqrt(FFcorfacerr.^2+var(nCl_arr))./bl;
    %%% save Cl
    flightCldat(ifield).l=l;
    flightCldat(ifield).bl=bl;
    flightCldat(ifield).Clflgiht=Clflight;
    flightCldat(ifield).nCl_arr=nCl_arr;
    flightCldat(ifield).FFcorfac=FFcorfac;
    flightCldat(ifield).FFcorfacerr=FFcorfacerr;
    flightCldat(ifield).meanmap=slopemean*cal;
    flightCldat(ifield).Cl_blcor=Cl_blcor;
    flightCldat(ifield).Cl_blcor_err=Cl_blcor_err;
end
save(strcat(savedir,'flightCldat'),'flightCldat');
return