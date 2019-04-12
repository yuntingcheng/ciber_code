function get_roughmask(flight,inst)
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

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

%%%%%%%%%%% get source mask %%%%%%%%%%%
% for ifield=4:8
%     fprintf('==================== ifild = %d ===================\n',ifield);
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

    %%%%%%% get maskin%%%%%%%
    [~,mask]=get_skymap(stackmapdat(ifield).cbmap,strmask.*mask_inst,4,5);

    %%%% some hand mask determined by eye(ex:Cosmic Rays) %%%%%%
    if inst==1 & ifield==8
    mask=circular_mask(442,40,40,mask);
    mask=circular_mask(742,342,10,mask);
    end

    if inst==2 & ifield==8
    mask=circular_mask(450,395,10,mask);
    mask=circular_mask(0,800,60,mask);
    end

    if inst==2 & ifield==5
    mask=elliptical_mask(215,220,60,15,80,mask);
    mask=elliptical_mask(710,90,80,20,60,mask);
    end
    
    maskdat.mask(ifield).bigmask=mask;
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