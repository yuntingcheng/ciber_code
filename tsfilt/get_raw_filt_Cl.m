function get_raw_filt_Cl(flight,inst)
pixscale=7;
mypaths=get_paths(flight);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/');
%%
for ifield=1:8

disp(sprintf('get Cl, ifield=%d',ifield));

dt=get_dark_times(flight,inst,ifield);
nfr_arr=2:dt.nfrhalf;
loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
load(strcat(loaddir,'bigmask'),'bigmask');
load(strcat(loaddir,'flightmap'),'flightmap');

%%%%%%%%%%%%%%%%% get Fourier Weight %%%%%%%%%%%%%%%%%
disp(sprintf('get Cl, ifield=%d -- get FW',ifield));

r2Cl2dstack_arr=zeros(numel(nfr_arr),1024,1024);
r1Cl2dstack_arr=zeros(numel(nfr_arr),1024,1024);
f2Cl2dstack_arr=zeros(numel(nfr_arr),1024,1024);
f1Cl2dstack_arr=zeros(numel(nfr_arr),1024,1024);

for i=1:numel(dt.time)
    load(strcat(loaddir,'labmap',num2str(i)),'labmap');
    for infr=1:numel(nfr_arr)
        rawmap=squeeze(labmap.rawmap_arr(infr,:,:));
        [~,maskin1]=get_skymap(rawmap,bigmask,4,5);
        rawmap=rawmap-mean(rawmap(find(maskin1)));
        rawmap=rawmap.*maskin1;
        rawmap=dc_offset_remove(rawmap,maskin1).*maskin1;
        [~,~,~,~,~,~,rCl2d] = get_angular_spec(rawmap,rawmap,pixscale);
        r2Cl2dstack_arr(infr,:,:)=squeeze(r2Cl2dstack_arr(infr,:,:))+rCl2d.^2;
        r1Cl2dstack_arr(infr,:,:)=squeeze(r1Cl2dstack_arr(infr,:,:))+rCl2d;

        
        filtmap=squeeze(labmap.filtmap_arr(infr,:,:));
        [~,maskin1]=get_skymap(filtmap,bigmask,4,5);
        filtmap=filtmap-mean(filtmap(find(maskin1)));
        filtmap=filtmap.*maskin1;
        filtmap=dc_offset_remove(filtmap,maskin1).*maskin1;
        [~,~,~,~,~,~,fCl2d] = get_angular_spec(filtmap,filtmap,pixscale);
        f2Cl2dstack_arr(infr,:,:)=squeeze(f2Cl2dstack_arr(infr,:,:))+fCl2d.^2;
        f1Cl2dstack_arr(infr,:,:)=squeeze(f1Cl2dstack_arr(infr,:,:))+fCl2d;
    end    
end

for infr=1:numel(nfr_arr)
    r1Cl2dstack=squeeze(r1Cl2dstack_arr(infr,:,:))./numel(dt.time);
    r2Cl2dstack=squeeze(r2Cl2dstack_arr(infr,:,:))./numel(dt.time);
    f1Cl2dstack=squeeze(f1Cl2dstack_arr(infr,:,:))./numel(dt.time);
    f2Cl2dstack=squeeze(f2Cl2dstack_arr(infr,:,:))./numel(dt.time);
    rCl2d_std=sqrt(r2Cl2dstack-r1Cl2dstack.^2);
    fCl2d_std=sqrt(f2Cl2dstack-f1Cl2dstack.^2);
    fw_raw=(fftshift(fftshift(1./squeeze(rCl2d_std))))';
    fw_filt=(fftshift(fftshift(1./squeeze(fCl2d_std))))';
    fwdat(infr).fw_raw=fw_raw;
    fwdat(infr).fw_filt=fw_filt;
end

clear r1Cl2dstack_arr r2Cl2dstack_arr f1Cl2dstack_arr f2Cl2dstack_arr

%%%%%%%%%%%%%%%%% get Cl of dark diff %%%%%%%%%%%%%%%%%
disp(sprintf('get Cl, ifield=%d -- get Cl dark',ifield));

rCl_arr=zeros(numel(dt.time),numel(nfr_arr),29);
fCl_arr=zeros(numel(dt.time),numel(nfr_arr),29);
wrCl_arr=zeros(numel(dt.time),numel(nfr_arr),29);
wfCl_arr=zeros(numel(dt.time),numel(nfr_arr),29);
for i=1:numel(dt.time)
    load(strcat(loaddir,'labmap',num2str(i)),'labmap');
    for infr=1:numel(nfr_arr)
        rawmap=squeeze(labmap.rawmap_arr(infr,:,:));
        [~,maskin1]=get_skymap(rawmap,bigmask,4,5);
        rawmap=rawmap-mean(rawmap(find(maskin1)));
        rawmap=rawmap.*maskin1;
        rawmap=dc_offset_remove(rawmap,maskin1).*maskin1;
        
        [rCl,l] = get_angular_spec(rawmap,rawmap,pixscale);
        rCl_arr(i,infr,:)=rCl;
        wrCl = get_angular_spec(rawmap,rawmap,pixscale,'w',fwdat(infr).fw_raw);
        wrCl_arr(i,infr,:)=wrCl;
        
        filtmap=squeeze(labmap.filtmap_arr(infr,:,:));
        [~,maskin1]=get_skymap(filtmap,bigmask,4,5);
        filtmap=filtmap-mean(filtmap(find(maskin1)));
        filtmap=filtmap.*maskin1;
        filtmap=dc_offset_remove(filtmap,maskin1).*maskin1;

        fCl = get_angular_spec(filtmap,filtmap,pixscale);
        fCl_arr(i,infr,:)=fCl;
        wfCl = get_angular_spec(filtmap,filtmap,pixscale,'w',fwdat(infr).fw_filt);
        wfCl_arr(i,infr,:)=wfCl;
    end    
end
diffCldat.l=l;
diffCldat.dark(ifield).nfr_arr=nfr_arr;
diffCldat.dark(ifield).rCl_arr=rCl_arr;
diffCldat.dark(ifield).fCl_arr=fCl_arr;
diffCldat.dark(ifield).wrCl_arr=wrCl_arr;
diffCldat.dark(ifield).wfCl_arr=wfCl_arr;

%%%%%%%%%%%%%%%%% get Cl flight diff %%%%%%%%%%%%%%%%%
disp(sprintf('get Cl, ifield=%d -- get Cl flight',ifield));

rCl_arr=zeros(numel(nfr_arr),29);
fCl_arr=zeros(numel(nfr_arr),29);
wrCl_arr=zeros(numel(nfr_arr),29);
wfCl_arr=zeros(numel(nfr_arr),29);

for infr=1:numel(nfr_arr)
    rawmap=squeeze(flightmap.rawmap_arr(infr,:,:));
    [~,maskin1]=get_skymap(rawmap,bigmask,4,5);
    rawmap=rawmap-mean(rawmap(find(maskin1)));
    rawmap=rawmap.*maskin1;
    rawmap=dc_offset_remove(rawmap,maskin1).*maskin1;
    
    [rCl,l] = get_angular_spec(rawmap,rawmap,pixscale);
    rCl_arr(infr,:)=rCl;
    wrCl = get_angular_spec(rawmap,rawmap,pixscale,'w',fwdat(infr).fw_raw);
    wrCl_arr(infr,:)=wrCl;

    filtmap=squeeze(flightmap.filtmap_arr(infr,:,:));
    [~,maskin1]=get_skymap(filtmap,bigmask,4,5);
    filtmap=filtmap-mean(filtmap(find(maskin1)));
    filtmap=filtmap.*maskin1;
    filtmap=dc_offset_remove(filtmap,maskin1).*maskin1;
        
    fCl = get_angular_spec(filtmap,filtmap,pixscale);
    fCl_arr(infr,:)=fCl;
    wfCl = get_angular_spec(filtmap,filtmap,pixscale,'w',fwdat(infr).fw_filt);
    wfCl_arr(infr,:)=wfCl;   
end
diffCldat.flight(ifield).nfr_arr=nfr_arr;
diffCldat.flight(ifield).rCl_arr=rCl_arr;
diffCldat.flight(ifield).fCl_arr=fCl_arr;
diffCldat.flight(ifield).wrCl_arr=wrCl_arr;
diffCldat.flight(ifield).wfCl_arr=wfCl_arr;
save(sprintf('%s/diffCldat',savedir),'diffCldat');
end

return 