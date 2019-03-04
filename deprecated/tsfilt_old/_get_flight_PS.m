%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the diff PS with best fit Chi2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1=-3.6;
load('/Volumes/HD1TB/CIBER/tsfilt/sinfiltamp/DarkLong/DCtemplate',...
                                                    'DCtemplate');

savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
savedirhalf='/Users/ytcheng/ciber/doc/20170209_TsFilter/halfPS/';

load(strcat(savedirfull,'FFdat'),'FFdat');

load(strcat(savedirfull,'fwfulldat'),'fwfulldat');
load(strcat(savedirhalf,'fwhalfdat'),'fwhalfdat');

simdir='/Users/ytcheng/ciber/doc/20170320_FFsim/';
load(strcat(simdir,'simCldat'),'simCldat');
%%
for ifield=[8,7,6,5,4,2,1]
disp(sprintf('ifield=%d',ifield));
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
load(strcat(loaddir,'maskin'),'maskin');
load(strcat(loaddir,'flightmap'),'flightmap');
dt=get_dark_times(flight,inst,ifield);

slopemap=flightmap.rawmapf;
slopemean=mean(slopemap(find(maskin)));

FF=FFdat(ifield).FF;
FFcorfac=simCldat.auto(ifield).meanCld;
FFcorfacerr=simCldat.auto(ifield).stdCld;
%% get weighted noise PS plot Cl-l

disp(sprintf('get weighted Cl'));

nClf_arr=zeros(numel(dt.time),21);
nCld_arr=zeros(numel(dt.time),21);
for i=1:numel(dt.time)
    i
    load(strcat(loaddir,'labdat',num2str(i)),'labdat');
    
    filtmapf=labdat.filtmapf;
    phmap=photonnoise_realization(slopemean.*ones(1024),G1,dt.nfr,frate);
    filtmapf=(filtmapf-DCtemplate)./FF+phmap; 
    W=(fftshift(fftshift(1./fwfulldat(ifield).std_fCl2d)))';
    [Clf]=get_Cl(filtmapf,maskin,fwfulldat(ifield).Mkk,pixscale,W);
    nClf_arr(i,:)=Clf.*cal.^2;
    
    filtmapd=(labdat.filtmap1-labdat.filtmap2)./FF./2;
    filtmapd(find(filtmapd~=filtmapd))=0;
    filtmapd(find(filtmapd==inf))=0;filtmapd(find(filtmapd==-inf))=0;
    maskd=maskin;maskd(find(filtmapd==0))=0;
    [~,maskd]=get_skymap(filtmapd,maskd,4,5);
    phmap1=photonnoise_realization(slopemean.*ones(1024),G1,dt.nfrhalf,frate);
    phmap2=photonnoise_realization(slopemean.*ones(1024),G1,dt.nfrhalf,frate);
    phmap=(phmap1-phmap2)./2;
    filtmapd=filtmapd+phmap;    
    W=(fftshift(fftshift(1./fwhalfdat(ifield).std_dCl2d)))';
    [Cld]=get_Cl(filtmapd,maskd,fwhalfdat(ifield).Mkk,pixscale,W);
    nCld_arr(i,:)=Cld.*cal.^2;
end
flightf=flightmap.filtmapf;
flightf=(flightf-DCtemplate)./FF.*maskin;
plane=plane_fit(flightf,maskin);
flightf=flightf-plane.*maskin;
flightf(find(flightf~=flightf))=0;
flightf(find(flightf==-inf))=0;
flightf(find(flightf==inf))=0;
calmapf=flightf.*cal;
flightmap.calmapf=calmapf;
W=(fftshift(fftshift(1./fwfulldat(ifield).std_fCl2d)))';
[Clflightf]=get_Cl(calmapf,maskin,fwfulldat(ifield).Mkk,pixscale,W);

flightd=(flightmap.filtmap1+flightmap.filtmap2)./FF./2;
flightd(find(flightd~=flightd))=0;
flightd(find(flightd==-inf))=0;
flightd(find(flightd==inf))=0;
calmapd=flightd.*cal;
flightmap.calmapd=calmapd;
W=(fftshift(fftshift(1./fwhalfdat(ifield).std_dCl2d)))';
[Clflightd,~,~,l]=get_Cl(calmapd,maskin,fwhalfdat(ifield).Mkk,pixscale,W);

Clflightf=Clflightf;
nClf_arr=nClf_arr;
Clflightd=Clflightd;
nCld_arr=nCld_arr;

save(strcat(loaddir,'flightmap'),'flightmap');
%% get bl 
load(strcat('/Volumes/HD1TB/CIBER/data/',num2str(flight),...
'/dr/TM',num2str(inst),'_',dt.name,'_dr150206.mat'));
logbl=interp1(log10(data.psf.l),log10(data.psf.bl),log10(l),'spline');
bl=10.^logbl;
%% save Cl
flightCldat(ifield).l=l;
flightCldat(ifield).bl=bl;
flightCldat(ifield).Clflgihtf=Clflightf;
flightCldat(ifield).nClf_arr=nClf_arr;
flightCldat(ifield).FFcorfac=FFcorfac;
flightCldat(ifield).FFcorfacerr=FFcorfacerr;
flightCldat(ifield).meanmap=slopemean*cal;
flightCldat(ifield).Cl_blcor=(Clflightf-FFcorfac-mean(nClf_arr))./bl;
flightCldat(ifield).Cl_blcor_err=sqrt(FFcorfacerr.^2+var(nClf_arr))./bl;
end
save(strcat(savedirfull,'flightCldat'),'flightCldat');
