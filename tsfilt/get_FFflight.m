function get_FFflight(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stack the flight FF
%The only field 4,5,6,7,8 is used. 
%The field itself is excluded. 
%ex: FF for field8 is stacking 4,5,6,7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',savedir),'maskdat');

for ifield=4:8
    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'flightmap'),'flightmap');

    mask=maskdat.mask(ifield).bigmask;
    meanmap=mean(flightmap.rawmapf(find(mask))).*cal;

    data(ifield).obs=flightmap.filtmapf.*cal;
    data(ifield).mask=mask;
    data(ifield).meanmap=meanmap;
end

for ifield=4:8
    FF=zeros(1024);stack_mask=zeros(1024);
    for jfield=4:8
        if jfield ~= ifield
            obs=data(jfield).obs;
            mask=data(jfield).mask;
            meanmap=data(jfield).meanmap;
            FF=FF+(obs.*mask./sqrt(meanmap));
            stack_mask=stack_mask+mask.*sqrt(meanmap);
        end
    end
FF=FF./stack_mask;FF((find(FF~=FF)))=0;
FFdat(ifield).FF=FF;

%%% make the unholy flat
FFmask=ones(1024);FFmask(find(FF==0))=0;
FFunholy = unholy_map(FF,FFmask,200);
FFdat(ifield).FFunholy = FFunholy;
end

save(strcat(savedir,'FFdat'),'FFdat');
return