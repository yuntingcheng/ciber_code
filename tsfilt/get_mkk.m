function get_mkk(flight,inst)

pixscale=7;
mypaths=get_paths(flight);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',savedir),'maskdat');
load(sprintf('%s/fwdat',savedir),'fwdat');

[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

for ifield=8:-1:4
    bigmask=maskdat.mask(ifield).bigmask;
    
    fw=fwdat(ifield).fw_filt;

    disp(sprintf('get Mkk, field%d weighted full',ifield));
    mkk_full =get_mkk_sim(bigmask,pixscale,binl,100,numel(binl),1,fw,0,NaN);
                                                      
    disp(sprintf('get Mkk, field%d unweighted',ifield));
    mkk_nw =get_mkk_sim(bigmask,pixscale,binl,100,numel(binl),1,ones(1024),0,NaN);
    
    mkkdat.auto(ifield).mask=bigmask;
    mkkdat.auto(ifield).mkk_wfull=mkk_full;
    mkkdat.auto(ifield).mkk_nw=mkk_nw;
save(sprintf('%s/mkkdat',savedir),'mkkdat');
end
save(sprintf('%s/mkkdat',savedir),'mkkdat');
return