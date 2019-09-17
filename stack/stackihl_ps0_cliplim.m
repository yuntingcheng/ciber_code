function [clipmax_arr, clipmin_arr,rbins]=...
stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
mask_inst,strnum,Nsrc,verbose,idx_stack_arr,rmin,use2m)
%%
if use2m
    srcdat = tm_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
    'sample_type','all');
    subx_arr = srcdat.xs_arr;
    suby_arr = srcdat.ys_arr;
    subm_arr = srcdat.ms_arr;

else
    srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
    'sample_type','all');
    subx_arr = [srcdat.xg_arr,srcdat.xs_arr];
    suby_arr = [srcdat.yg_arr,srcdat.ys_arr];
    subm_arr = [srcdat.mg_arr,srcdat.ms_arr];

end

mask_inst = squeeze(mask_inst(inst,:,:));

%%% set up stacking %%%
if isnan(rmin)
    rad_arr = get_mask_radius_th(inst,ifield,subm_arr,0.5);
else
    rad_arr = get_mask_radius_th(inst,ifield,subm_arr,0.5,'rmin',rmin);
end

idx_arr = 1:numel(subm_arr);

if Nsrc ~= 0 & Nsrc < numel(idx_arr)
    idx_arr = datasample(idx_arr,Nsrc,'Replace',false);
end

if numel(idx_arr)>20
    idx_print=floor(numel(idx_arr)/20) * (1:20);
else
    idx_print=[];
end
print_count = 0;

if numel(idx_stack_arr) == 0
    idx_stack_arr = idx_arr;
end

nbins = 25;
dx = 1200;
profile = radial_prof(ones(2*dx+1),ones(2*dx+1),dx+1,dx+1,1,nbins);
rbinedges = profile.binedges;
rbins = binedges2bins(rbinedges).*0.7;

for ibin=1:nbins
    dat(ibin).cb = [];
    dat(ibin).ps = [];
end
ii=0;
for i=idx_stack_arr
    ii=ii+1;
    %%% print
    if ismember(ii,idx_print) & verbose
        print_count = print_count + 1;
        disp(sprintf('find cliplim %d %% %d/%d src between %.1f<m<%.1f '...
            ,print_count*5,ii,numel(idx_stack_arr),m_min,m_max));         
    end

    %%% unmask the target
    radmap = make_radius_map(cbmap,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
    if numel(sp1)==0
        continue
    end
    ri = radmap(sp1).*10;%%% subpixel unit
    cbi = cbmap(sp1);
    psi = psmap(sp1);
    inbins = find(rbinedges > max(ri));
    inbins = inbins(1)-1;
    
    for ibin = 1:inbins
        sp = find((ri>=rbinedges(ibin)) & (ri<rbinedges(ibin+1)));
        dat(ibin).cb = [dat(ibin).cb, cbi(sp)'];
        dat(ibin).ps = [dat(ibin).ps, psi(sp)'];
    end
    
end
%%
clipmin_arr = ones([2,nbins]).*-inf;
clipmax_arr = ones([2,nbins]).*inf;


Q1 = quantile([dat(1).cb,dat(2).cb,dat(3).cb,dat(4).cb],0.25);
Q3 = quantile([dat(1).cb,dat(2).cb,dat(3).cb,dat(4).cb],0.75);
IQR = Q3 - Q1;
for ibin=1:4
    clipmax_arr(1,ibin) = Q3+3*IQR;
    clipmin_arr(1,ibin) = Q1-3*IQR;
end
Q1 = quantile([dat(1).ps,dat(2).ps,dat(3).ps,dat(4).ps],0.25);
Q3 = quantile([dat(1).ps,dat(2).ps,dat(3).ps,dat(4).ps],0.75);
IQR = Q3 - Q1;
for ibin=1:4
    clipmax_arr(2,ibin) = Q3+3*IQR;
    clipmin_arr(2,ibin) = Q1-3*IQR;
end


for ibin=5:nbins
    if numel(dat(ibin).cb) == 0
        continue
    end
    
    Q1 = quantile(dat(ibin).cb,0.25);
    Q3 = quantile(dat(ibin).cb,0.75);
    IQR = Q3 - Q1;
    clipmax_arr(1,ibin) = Q3+3*IQR;
    clipmin_arr(1,ibin) = Q1-3*IQR;
    
    Q1 = quantile(dat(ibin).ps,0.25);
    Q3 = quantile(dat(ibin).ps,0.75);
    IQR = Q3 - Q1;
    clipmax_arr(2,ibin) = Q3+3*IQR;
    clipmin_arr(2,ibin) = Q1-3*IQR;

end

return
