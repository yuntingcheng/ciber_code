%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the diff PS with best fit Chi2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;
frate=cp(inst).framerate;
G1=-3.6;
load('/Volumes/HD1TB/CIBER/tsfilt/sinfiltamp/DarkLong/DCtemplate',...
                                                    'DCtemplate');

savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
savedirhalf='/Users/ytcheng/ciber/doc/20170209_TsFilter/halfPS/';
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
%% get Fourier weight

disp(sprintf('get Fourier weight'));
fCl2d_arr=zeros(numel(dt.time),1024,1024);
dCl2d_arr=zeros(numel(dt.time),1024,1024);
fCl2dph_arr=zeros(numel(dt.time),1024,1024);
dCl2dph_arr=zeros(numel(dt.time),1024,1024);

for i=1:numel(dt.time)
    load(strcat(loaddir,'labdat',num2str(i)),'labdat');
    
    filtmapf=labdat.filtmapf;
    filtmapf=(filtmapf-DCtemplate).*maskin;
    filtmapf=filtmapf-mean(filtmapf(find(maskin)));
    filtmapf=filtmapf.*maskin;
    [~,~,~,~,binl,~,fCl2d] = get_angular_spec(filtmapf,filtmapf,pixscale);
    fCl2d_arr(i,:,:)=fCl2d;
    
    filtmapf=labdat.filtmapf;
    phmap=photonnoise_realization(slopemean.*ones(1024),G1,dt.nfr,frate);
    filtmapf=(filtmapf-DCtemplate+phmap).*maskin;    
    filtmapf=filtmapf-mean(filtmapf(find(maskin)));
    filtmapf=filtmapf.*maskin;
    [~,~,~,~,~,~,fCl2d] = get_angular_spec(filtmapf,filtmapf,pixscale);
    fCl2dph_arr(i,:,:)=fCl2d;
    
    filtmapd=(labdat.filtmap1-labdat.filtmap2)./2;
    [~,maskd]=get_skymap(filtmapd,maskin,4,5);
    filtmapd=filtmapd-mean(filtmapd(find(maskd)));
    filtmapd=filtmapd.*maskd;
    [~,~,~,~,~,~,dCl2d] = get_angular_spec(filtmapd,filtmapd,pixscale);
    dCl2d_arr(i,:,:)=dCl2d;

    filtmapd=(labdat.filtmap1-labdat.filtmap2)./2;
    [~,maskd]=get_skymap(filtmapd,maskin,4,5);
    phmap1=photonnoise_realization(slopemean.*ones(1024),G1,dt.nfrhalf,frate);
    phmap2=photonnoise_realization(slopemean.*ones(1024),G1,dt.nfrhalf,frate);
    phmap=(phmap1-phmap2)./2;
    filtmapd=(filtmapd+phmap).*maskd;    
    filtmapd=filtmapd-mean(filtmapd(find(maskd)));
    filtmapd=filtmapd.*maskd;
    [~,~,~,~,~,~,dCl2d] = get_angular_spec(filtmapd,filtmapd,pixscale);
    dCl2dph_arr(i,:,:)=dCl2d;
    
end
Wfrn=(fftshift(fftshift(1./squeeze(std(fCl2d_arr)))))';
Wdrn=(fftshift(fftshift(1./squeeze(std(dCl2d_arr)))))';
Wfph=(fftshift(fftshift(1./squeeze(std(fCl2dph_arr)))))';
Wdph=(fftshift(fftshift(1./squeeze(std(dCl2dph_arr)))))';
fwfulldat(ifield).avg_fCl2d=squeeze(mean(fCl2d_arr));
fwfulldat(ifield).std_fCl2d=squeeze(std(fCl2d_arr));
fwhalfdat(ifield).avg_dCl2d=squeeze(mean(dCl2d_arr));
fwhalfdat(ifield).std_dCl2d=squeeze(std(dCl2d_arr));

%% get Mkk
disp(sprintf('get Mkk'));

fwfulldat(ifield).mask=maskin;
Mkk =get_mkk_sim(maskin,pixscale,binl,100,numel(binl),1,Wfrn,0,NaN);
fwfulldat(ifield).Mkk=Mkk;

fwhalfdat(ifield).mask=maskin;
Mkk =get_mkk_sim(maskin,pixscale,binl,100,numel(binl),1,Wdrn,0,NaN);
fwhalfdat(ifield).Mkk=Mkk;
end

save(strcat(savedirfull,'fwfulldat'),'fwfulldat');
save(strcat(savedirhalf,'fwhalfdat'),'fwhalfdat');
