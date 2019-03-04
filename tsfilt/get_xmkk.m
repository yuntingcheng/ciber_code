function get_xmkk(flight,inst)
mypaths=get_paths(flight);
pixscale=7;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',savedir),'maskdat');
load(sprintf('%s/fwdat',savedir),'fwdat');
load(sprintf('%s/mkkdat',savedir),'mkkdat');

[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

count=0;

for ifield=8:-1:1
    for jfield=8:-1:1
        if ifield>jfield
            count=count+1;
            disp(sprintf('count=%d,[i,j]=[%d,%d]',count,ifield,jfield));
            
            mask=maskdat.mask(ifield).bigmask.*maskdat.mask(jfield).bigmask;
            fw=(fwdat(ifield).fw_filt+fwdat(jfield).fw_filt)./2;
            mkk=get_mkk_sim(mask,pixscale,binl,100,numel(binl),1,fw,0,NaN);
                            
            mkkdat.cross(count).field=[ifield jfield];
            mkkdat.cross(count).mask=mask;
            mkkdat.cross(count).mkk=mkk;
            save(sprintf('%s/mkkdat',savedir),'mkkdat');
        end
    end
end
save(sprintf('%s/mkkdat',savedir),'mkkdat');
return