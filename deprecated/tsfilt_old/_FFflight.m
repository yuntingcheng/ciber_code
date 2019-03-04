%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stack FF flight
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
%% load in the infomation
for ifield=[8,7,6,5,4,2,1]
    ifield
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
load(strcat(loaddir,'maskin'),'maskin');
load(strcat(loaddir,'flightmap'),'flightmap');

meanmap=mean(flightmap.rawmapf(find(maskin))).*cal;

data(ifield).obs=flightmap.filtmapf.*cal;
data(ifield).mask=maskin;
data(ifield).meanmap=meanmap;
end
%%
for ifield=[8,7,6,5,4,2,1]
FF=zeros(1024);stack_mask=zeros(1024);
FFuse=[0,0,0,1,1,1,1,1];FFuse(ifield)=0;
for jfield=[8,7,6,5,4,2,1]
    if FFuse(jfield)==1
    obs=data(jfield).obs;
    mask=data(jfield).mask;
    meanmap=data(jfield).meanmap;
    FF=FF+(obs.*mask./sqrt(meanmap));
    stack_mask=stack_mask+mask.*sqrt(meanmap);
    end
end
FF=FF./stack_mask;FF((find(FF~=FF)))=0;
FFdat(ifield).FF=FF;
end
%%
savedir='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
save(strcat(savedir,'FFdat'),'FFdat');