function get_fweight(flight,inst)
mypaths=get_paths(flight);
pixscale=7;

%%% get DC template
load(sprintf('%sTM%d/darklongdat',mypaths.filtmap,inst),'darklongdat');
DCtemplate=darklongdat.DCtemplate;

maskdir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',maskdir),'maskdat');

for ifield=4:8
    disp(sprintf('get fweight, ifield=%d',ifield));
    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    bigmask=maskdat.mask(ifield).bigmask;
    dt=get_dark_times(flight,inst,ifield);

    rCl2d_arr=zeros(numel(dt.time),1024,1024);
    fCl2d_arr=zeros(numel(dt.time),1024,1024);

    for i=1:numel(dt.time)
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');

        %%% fw of unfiltered full %%%%
        rawmapf=labmap.rawmapf;
        rawmapf=(rawmapf-DCtemplate).*bigmask;
        rawmapf=rawmapf-mean(rawmapf(bigmask==1));
        bigmask1 = sigclip_mask(rawmapf,bigmask,5,3);
        rawmapf=rawmapf.*bigmask1;
        rawmapf=dc_offset_remove(rawmapf,bigmask1).*bigmask1;
        [~,~,~,~,~,~,rCl2d] = get_angular_spec(rawmapf,rawmapf,pixscale);
        rCl2d_arr(i,:,:)=rCl2d;

        %%% fw of filtered full %%%%
        filtmapf=labmap.filtmapf;
        filtmapf=(filtmapf-DCtemplate).*bigmask;
        filtmapf=filtmapf-mean(filtmapf(bigmask==1));
        bigmask1 = sigclip_mask(filtmapf,bigmask,5,3);
        filtmapf=filtmapf.*bigmask1;
        filtmapf=dc_offset_remove(filtmapf,bigmask1).*bigmask;
        [~,~,~,~,~,~,fCl2d] = get_angular_spec(filtmapf,filtmapf,pixscale);
        fCl2d_arr(i,:,:)=fCl2d;

    end

    fw_raw=(fftshift(fftshift(1./squeeze(std(rCl2d_arr)))))';
    fw_filt=(fftshift(fftshift(1./squeeze(std(fCl2d_arr)))))';
    fwdat(ifield).fw_raw=fw_raw;
    fwdat(ifield).fw_filt=fw_filt;
end
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%sfwdat',savedir),'fwdat');
return