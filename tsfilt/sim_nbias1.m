function sim_nbias1(flight,inst,G1)
%%
pixscale=7;
mypaths=get_paths(flight);
cp=get_cal_params('flight',flight);
frate=cp(inst).framerate;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

load(strcat(savedir,'darkstat'),'darkstat');
load(strcat(savedir,'fwdat'),'fwdat');
load(strcat(savedir,'mkkdat'),'mkkdat');

%%% load labflat                
dataname = strcat(mypaths.release,'TM',num2str(inst),'_SWIRE_dr150206.mat');
% lab flat are the same for all fields
load(dataname);labflat=data.flat.labflat;
clear data
%% simulation input info
for ifield=4:8
    siminfo(ifield).weight=fwdat(ifield).fw_filt;
    siminfo(ifield).Cl2d_ave=darkstat(ifield).Cl2d_std;
    siminfo(ifield).bigmask=mkkdat.auto(ifield).mask;
    siminfo(ifield).Mkk=mkkdat.auto(ifield).mkk_wfull;
end
clear darkstat fwdat mkkdat
%%
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);

    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'flightmap'),'flightmap');
    slopemap=flightmap.rawmapf;
    meanmap=mean(slopemap(find(siminfo(ifield).bigmask)));
    clear flightfr slopemap

    siminfo(ifield).name=dt.name;
    siminfo(ifield).nfr=dt.nfr;
    siminfo(ifield).meanmap=meanmap;
end
%% 100 runs
tic
nsim=100;
nCl_arr=zeros(8,nsim,21);
% rnClx_arr=zeros(size(mkkdat.cross,2),nsim,21);
% phClx_arr=zeros(size(mkkdat.cross,2),nsim,21);
%%
for isim=1:nsim
%%
%%%%% get sim maps %%%%%%
for ifield =4:8
    meanmap=siminfo(ifield).meanmap;
    bigmask=siminfo(ifield).bigmask;
    slopemap=(meanmap).*labflat;
    
    phmap=photonnoise_realization(slopemap,G1,siminfo(ifield).nfr,frate);
    rnmap=readnoise_realization(siminfo(ifield).Cl2d_ave,pixscale,'norand',1);
    
    nmap=rnmap + phmap;
    nobs=(slopemap+phmap+rnmap).*labflat.*bigmask;
    
    mapdat(ifield).nmap=nmap;
    mapdat(ifield).nobs=nobs;
end

%%%% stack FF %%%%%
for ifield=4:8
nFF=zeros(size(labflat));
nstack_mask=zeros(size(labflat));

FFuse=[0,0,0,1,1,1,1,1];FFuse(ifield)=0;
for jfield=4:8
    if FFuse(jfield)==1
        bigmask=siminfo(jfield).bigmask;

        nobs=mapdat(jfield).nobs;
        mean_n=mean(nobs(find(nobs)));
        nFF=nFF+(nobs./sqrt(mean_n));
        nstack_mask=nstack_mask+bigmask.*sqrt(mean_n);
    end
end
nFF=nFF./nstack_mask;nFF((find(nFF~=nFF)))=0;
mapdat(ifield).nFF=nFF;
end

%%%%% get FF corrected map %%%%%
for ifield=4:8
    mapdat(ifield).fnmap=mapdat(ifield).nmap;
    nmap1 = mapdat(ifield).nobs./mapdat(ifield).nFF;
    bigmask1 = sigclip_mask(nmap1,bigmask,5,5);
    mapdat(ifield).nmap1=mapdat(ifield).nobs./mapdat(ifield).nFF;
    mapdat(ifield).bigmask1 = bigmask1;
end

%%%%%% get auto spectra %%%%%%
for ifield=4:8
    bigmask=mapdat(ifield).bigmask1;
    fw=siminfo(ifield).weight;
    Mkk=siminfo(ifield).Mkk;    
    
    [nCl]=get_Cl(mapdat(ifield).fnmap,bigmask,Mkk,pixscale,fw);
    [nCl1,~,~,l]=get_Cl(mapdat(ifield).nmap1,bigmask,Mkk,pixscale,fw);
    
    nCl_arr(ifield,isim,:)=nCl1-nCl;
end

cal = -170.3608;
loglog(l,l.*(l+1).*nCl.*cal.*cal./2./pi,'k');hold on
loglog(l,l.*(l+1).*nCl1.*cal.*cal./2./pi,'r');hold on

drawnow
%%
%%%%%% get cross spectra %%%%%%
% count=0;
% for ifield=8:-1:4
%     bigmaski=siminfo(ifield).bigmask;
%     fwi=siminfo(ifield).weight;
%     rnmap1i=mapdat(ifield).rnmap1;
%     phmap1i=mapdat(ifield).phmap1;
% 
%     for jfield=8:-1:4
%         if jfield<ifield
%             count=count+1;
%             bigmaskj=siminfo(jfield).bigmask;
%             bigmask=bigmaski.*bigmaskj;
%             fwj=siminfo(jfield).weight;
%             fw=(fwi+fwj)./2;
%             Mkk=mkkdat.cross(count).mkk;
%             rnmap1j=mapdat(jfield).rnmap1;
%             phmap1j=mapdat(jfield).phmap1;
%             [rnClx]=get_Clx(rnmap1i,rnmap1j,bigmask,Mkk,pixscale,fw);
%             [phClx]=get_Clx(phmap1i,phmap1j,bigmask,Mkk,pixscale,fw);
%             rnClx_arr(count,isim,:)=rnClx;
%             phClx_arr(count,isim,:)=phClx;
%         end        
%     end   
% end

%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('isim=%d,t=%.1f min',isim,toc/60));
%%
end
%%
simCldat.l=l;
for ifield=4:8
    navgCld=squeeze(mean(nCl_arr(ifield,:,:)));
    nstdCld=squeeze(std(nCl_arr(ifield,:,:)));

    simCldat.nbias.auto(ifield).navgCl=navgCld';
    simCldat.nbias.auto(ifield).nstdCl=nstdCld';
end
%% cross spectrum
% count=0;
% for ifield=8:-1:4
%     for jfield=8:-1:4
%         if jfield<ifield
%             count=count+1;
%             rnavgClx=squeeze(mean(rnClx_arr(count,:,:)));
%             rnstdClx=squeeze(std(rnClx_arr(count,:,:)));
%             phavgClx=squeeze(mean(phClx_arr(count,:,:)));
%             phstdClx=squeeze(std(phClx_arr(count,:,:)));
%             
%             simCldat.nbias.cross(count).field=mkkdat.cross(count).field;
%             simCldat.nbias.cross(count).rnavgClx=rnavgClx';
%             simCldat.nbias.cross(count).rnstdClx=rnstdClx';
%             simCldat.nbias.cross(count).phavgClx=phavgClx';
%             simCldat.nbias.cross(count).phstdClx=phstdClx';
%         end
%     end
% end
%%
save(strcat(savedir,'simCldat'),'simCldat');
return