function get_bigmask(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get bigmask which include:
% 1. maskin: mask used in time stream filter, 
%            including MZ mask and crazy pix mask for all darks.
% 2. jackmask:3 sigma clip of two halves diff
% 3. hand mask: mask crazy pix in diff map by hands,
%               these are cosmic rays or residual PSF.
% 4. FFmask: pix don't have FF. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

maskinstdir=strcat(mypaths.ciberdir,'doc/20160808_DarkProcess/40030/');
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',savedir),'maskdat');

%%%%%%%%%%% get source mask %%%%%%%%%%%
% for ifield=4:8
%     fprintf('==================== ifild = %d ===================\n',ifield);
%     %%% MZ mask
%     %load(sprintf('%sTM%d_%s_dr150206.mat',mypaths.release,inst,dt.name));
%     %maskmz=double(~data.mask.mask);
% 
%     [tmmask,tmnum] = make_mask_2m(flight,inst,'y',ifield,1,13);
%     [psmask,psnum] = make_mask_ps(flight,inst,'y',ifield,0,0,23);
%     strmask = psmask.*tmmask;
%     strnum = psnum + tmnum;
%     
%     maskdat.mask(ifield).strmask=strmask;
%     maskdat.mask(ifield).strnum=strnum;
% end
% save(sprintf('%s/maskdat',savedir),'maskdat');


load(sprintf('%s/maskdat',savedir),'maskdat');
for ifield=4:8    
    %%% mask_inst
    maskdat.mask(ifield).mask_inst=mask_inst;
    
    strmask = maskdat.mask(ifield).strmask;

    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'flightmap'),'flightmap');

    %%%%%%% get maskin%%%%%%%
    [~,maskin]=get_skymap(flightmap.rawmapf,strmask.*mask_inst,4,5);
    [~,maskin]=get_skymap(flightmap.filtmapf,maskin,4,5);

    %%%%%%% jack mask %%%%%%
    diffmap=flightmap.rawmap1-flightmap.rawmap2;
    [~,jackmask]=get_skymap(diffmap,maskin,3,5);

    %%%% some hand mask determined by eye(ex:Cosmic Rays) %%%%%%
    if inst==1 & ifield==8
    jackmask=circular_mask(442,40,40,jackmask);
    jackmask=circular_mask(742,342,10,jackmask);
    end

    if inst==2 & ifield==8
    jackmask=circular_mask(450,395,10,jackmask);
    jackmask=circular_mask(0,800,60,jackmask);
    end

    if inst==2 & ifield==5
    jackmask=elliptical_mask(215,220,60,15,80,jackmask);
    jackmask=elliptical_mask(710,90,80,20,60,jackmask);
    end
    
    masktemp(ifield).mask=jackmask;
end

%%%%%%%%%%% get FF mask %%%%%%%%%%%
for ifield=4:8
    FFmask=zeros(1024);
    for jfield=4:8
        if jfield~=ifield
            FFmask=FFmask+masktemp(jfield).mask;
        end
        FFmask(find(FFmask))=1;
    end
    bigmask=FFmask.*masktemp(ifield).mask;
    savedir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    save(strcat(savedir,'bigmask'),'bigmask');

    maskdat.mask(ifield).bigmask=bigmask;
end

%%%%%%%%%%% exclude strmask for customize star mask %%%%
for ifield=4:8
    bigmask=maskdat.mask(ifield).bigmask;
    strmask=maskdat.mask(ifield).strmask;
    mask_inst=maskdat.mask(ifield).mask_inst;
    mask1=mask_inst.*strmask;
    sp = find(bigmask==0 & mask1==1);
    maskcrap=ones(1024);
    maskcrap(sp)=0;
    maskdat.mask(ifield).maskcrap=maskcrap;
    maskdat.mask(ifield).nosrc=maskcrap.*mask_inst;
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/maskdat',savedir),'maskdat');

return