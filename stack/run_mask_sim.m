function run_mask_sim(flight,inst,field,hsc_idx)

mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');

if strcmp(field,'SIDES')
    for ifield=4:8
        [mask,num]=make_mask_sim(flight,inst,ifield,-5,20);
        maskdat.mask(ifield).strmask_sim = mask;
        maskdat.mask(ifield).strnum_sim = num;
    end
elseif strcmp(field,'HSC')
    fieldname = HSC_field_name(hsc_idx);
    for ifield=4:8
        [mask,num]=make_mask_hsc(flight,inst,ifield,fieldname,-5,20);
        maskdat.mask(ifield).hsc(hsc_idx).strmask_sim = mask;
        maskdat.mask(ifield).hsc(hsc_idx).strnum_sim = num;
    end
end

save(strcat(loaddir,'maskdat'),'maskdat');
return