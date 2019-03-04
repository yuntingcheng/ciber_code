function run_mask_sim(flight,inst)

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');

for ifield=4:8
    [mask,num]=make_mask_sim(flight,inst,ifield,-5,20);
    maskdat.mask(ifield).strmask_sim = mask;
    maskdat.mask(ifield).strnum_sim = num;
end

save(strcat(loaddir,'maskdat'),'maskdat');
return