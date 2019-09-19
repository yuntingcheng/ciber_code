function [hscprof,hscerr] = get_hsc_clus_prof(flight, inst, masklim)

mypaths=get_paths(flight);

load(sprintf('%s/TM%d/psfdat',mypaths.alldat,inst),'psfdatallfields');
psf_arr = psfdatallfields(8).psfps;
r_arr = psfdatallfields(8).r_arr;
hscprof = zeros([4,numel(psf_arr)]);
hscerr = zeros([4,numel(psf_arr)]);
for im=1:4
    excb_all = zeros([12,numel(r_arr)]);
    exps_all = zeros([12,numel(r_arr)]);
    excb_err_all = zeros([12,numel(r_arr)]);
    exps_err_all = zeros([12,numel(r_arr)]);
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
    hscprof(im,:) = excb_all;
    hscerr(im,:) = excb_err_all;    
end

return