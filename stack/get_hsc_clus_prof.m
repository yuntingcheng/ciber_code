function hscclusdat = get_hsc_clus_prof(flight, inst, masklim, weight)

mypaths=get_paths(flight);

load(sprintf('%s/TM%d/psfdat_SWIRE',mypaths.alldat,inst),'psfdatall');
r_arr = psfdatall.comb(1).r_arr;
for im=1:4
    psf_arr = psfdatall.comb(im).all.profps;
    excb_all = zeros([12,numel(r_arr)]);
    excb_err_all = zeros([12,numel(r_arr)]);
    for hsc_idx=0:11
        [name,~] = HSC_fields_info(hsc_idx);
        loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
        if masklim
            load(sprintf('%s/hsc/stackdathsc_%s_masklim',...
                loaddir,name),'stackdathsc');  
            load(sprintf('%s/hsc/stackdathsc_%s_masklim_bk',...
                loaddir,name),'stackdathscbk');  
        else
            load(sprintf('%s/hsc/stackdathsc_%s',...
                loaddir,name),'stackdathsc');
            load(sprintf('%s/hsc/stackdathsc_%s_bk',...
                loaddir,name),'stackdathscbk');
        end
        
        stackdat = stackdathsc(im).stackdat;
        bkallcb_arr = ones([numel(stackdathscbk),numel(r_arr)]).*nan;
        for iter=1:numel(stackdathscbk)
            stackdatbk = stackdathscbk(iter).stackdat(im);
            sp  = find(stackdatbk.hitmapg_arr~=0);
            bkallcb_arr(iter,sp) = stackdatbk.profcbg_arr(sp);
        end
        bkavgcb_arr = nanmean(bkallcb_arr);
        sp = find(bkavgcb_arr==bkavgcb_arr);
        bkavgcb_arr = spline(r_arr(sp),bkavgcb_arr(sp),r_arr);
        bkavgcb_err = nanstd(bkallcb_arr);
        sp = find(bkavgcb_err==bkavgcb_err);
        bkavgcb_err = spline(r_arr(sp),bkavgcb_err(sp),r_arr);

        profcbg = stackdat.all.profcbg - bkavgcb_arr;
        profcbg_err = sqrt(stackdat.errjack.profcbg.^2 + bkavgcb_err.^2);
        excb_all(hsc_idx+1,:) = profcbg - psf_arr.*profcbg(1);
        excb_err_all(hsc_idx+1,:) = profcbg_err;
    end
    excb_all = sum(excb_all./(excb_err_all.^2))./sum(1./(excb_err_all.^2));
    excb_err_all = sqrt(1./sum(1./(excb_err_all.^2)));
    hscclusdat(im).r_arr = r_arr;
    hscclusdat(im).prof = excb_all;
    hscclusdat(im).err = excb_err_all;  
    
    %%% fit a line %%%
    linefitpar = polyfit(log10(r_arr(r_arr > 7)),log10(excb_all(r_arr > 7)),1);
    hscclusdat(im).linefitpar = linefitpar;
    model = 10.^polyval(linefitpar,log10(r_arr));
    hscclusdat(im).linefitmodel = model;
    
    %%% binning ot 15 ard bins %%%
    [rsub_arr,~] = profile_radial_binning(r_arr,weight,find(r_arr>100));
    hscclusdat(im).rsub_arr = rsub_arr;
    [prof15,~] = profile_radial_binning(excb_all,weight,find(r_arr>100));
    [var15,~] = profile_radial_binning(excb_err_all.^2,weight,find(r_arr>100));
    err15 = sqrt(var15);
    hscclusdat(im).rsub_arr = rsub_arr;
    hscclusdat(im).profsub = prof15;
    hscclusdat(im).errsub = err15;
    [model15,~] = profile_radial_binning(model,weight,find(r_arr>100));
    hscclusdat(im).linefitmodelsub = model15;
    
end

return