function [profcb_arr,profps_arr,hitmap_arr]=...
stackihl_ps0_hist_map_bk(dx,cbmap,psmap,mask_inst,strmask,N_arr,verbose,subpix)
%%
Npix=size(cbmap,1);
x_arr=Npix*rand(1,max(N_arr))+0.5;
y_arr=Npix*rand(1,max(N_arr))+0.5;

idx_arr = 1:numel(x_arr);

if numel(idx_arr)>20
    idx_print=floor(numel(idx_arr)/20) * (1:20);
else
    idx_print=[];
end
print_count = 0;

nbins = 25;
profcb_arr = zeros([1,nbins]);
profps_arr = zeros([1,nbins]);
hitmap_arr = zeros([1,nbins]);

profcb_mat = zeros([numel(N_arr),nbins]);
profps_mat = zeros([numel(N_arr),nbins]);
hitmap_mat = zeros([numel(N_arr),nbins]);

if subpix
    rad = make_radius_map(zeros(2*dx+1),dx+1,dx+1);
    profile = radial_prof(rad,ones(2*dx+1),dx+1,dx+1,1,nbins);
    rbinedges = profile.binedges;
else
    rad = make_radius_map(zeros(2*dx+1),dx+1,dx+1);
    profile = radial_prof(rad,ones(2*dx+1),dx+1,dx+1,1,nbins);
    rbinedges = profile.binedges./10;
    rad = make_radius_map(zeros(2*dx/10+1),dx/10+1,dx/10+1);
end

for i=idx_arr

    cbmapi = cbmap.*strmask.*mask_inst;
    psmapi = psmap.*strmask.*mask_inst;
    
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');  
    if subpix
        %%% rebin to 10x finer map
        stampcb0 = imresize(cbmapi,10,'nearest');
        stampps0 = imresize(psmapi,10,'nearest');
        %%% get the stamp
        xcent = round(x_arr(i)*10-4.5) + dx;
        ycent = round(y_arr(i)*10-4.5) + dx;
        stampcb = stampcb0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
        stampps = stampps0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    else
        %%% get the stamp
        xcent = round(x_arr(i)) + dx/10;
        ycent = round(y_arr(i)) + dx/10;
        stampcb = cbmapi(xcent-dx/10:xcent+dx/10,ycent-dx/10:ycent+dx/10);
        stampps = psmapi(xcent-dx/10:xcent+dx/10,ycent-dx/10:ycent+dx/10);        
    end

    maskstamp = zeros(size(stampcb));
    maskstamp(find(stampcb~=0)) = 1;
    for ibin = 1:nbins
        sp = find((maskstamp~=0) & (rad>=rbinedges(ibin)) ...
            & (rad<rbinedges(ibin+1)));
        profcb_arr(ibin) = profcb_arr(ibin) + sum(stampcb(sp));
        profps_arr(ibin) = profps_arr(ibin) + sum(stampps(sp));
        hitmap_arr(ibin) = hitmap_arr(ibin) + numel(sp);
    end
    
    if ismember(i,N_arr)
        Nidx = find(N_arr==i);
        profcb_mat(Nidx,:) = profcb_arr./hitmap_arr;
        profps_mat(Nidx,:) = profps_arr./hitmap_arr;
        hitmap_mat(Nidx,:) = hitmap_arr;
    end

    %%% print
    if ismember(i,idx_print) & verbose
        print_count = print_count + 1;
        disp(sprintf('stack %d %% %d/%d src',...
            print_count*5,i,numel(idx_arr)));         
    end
end

profcb_arr = profcb_mat;
profps_arr = profps_mat;
hitmap_arr = hitmap_mat;

return