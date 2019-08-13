function [histbkdat]=stackihl_ps_randomN_hist_subreg...
    (N_arr,dx,cbmap,psmap,Ibinedges_cb,Ibinedges_ps,mask,isubreg,spire)

% randomize the position
sp = [];
while numel(sp)<max(N_arr)
    subx_arr=1023*rand(1,max(N_arr)*100)+1;
    suby_arr=1023*rand(1,max(N_arr)*100)+1;
    
    regidx = zeros(size(subx_arr));
    regidx((subx_arr<=512.5) & (suby_arr<=512.5)) = 1;
    regidx((subx_arr<=512.5) & (suby_arr>512.5)) = 3;
    regidx((subx_arr>512.5) & (suby_arr<=512.5)) = 5;
    regidx((subx_arr>512.5) & (suby_arr>512.5)) = 7;
    regidx((regidx==1) & (subx_arr>120.5) & (suby_arr>120.5)) = 2;
    regidx((regidx==3) & (subx_arr>120.5) & (suby_arr<=904.5)) = 4;
    regidx((regidx==5) & (subx_arr<=904.5) & (suby_arr>120.5)) = 6;
    regidx((regidx==7) & (subx_arr<=904.5) & (suby_arr<=904.5)) = 8;
    
    sp=find(regidx==isubreg);
    subx_arr=subx_arr(sp(1:max(N_arr)));
    suby_arr=suby_arr(sp(1:max(N_arr)));
end

nbins = 25;
histcb = zeros(nbins,numel(Ibinedges_cb)-1);
histps = zeros(nbins,numel(Ibinedges_ps)-1);

rad = make_radius_map(zeros(2*dx+1),dx+1,dx+1);
profile = radial_prof(rad,ones(2*dx+1),dx+1,dx+1,1,nbins);
binedges = profile.binedges;

idx_print=floor(max(N_arr)/100) * (1:100);
print_count = 0;

for i=1:max(N_arr)
    cbmapi = cbmap.*mask;
    psmapi = psmap.*mask;
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');
    %%% rebin to 10x finer map
    stampcb = imresize(cbmapi,10,'nearest');
    stampps = imresize(psmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stampcb = stampcb(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stampps = stampps(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp=zeros(size(stampcb));
    maskstamp(find(stampcb~=0))=1;
    
    for ibin = 1:nbins
        sp = find((stampcb~=0) & (rad>=binedges(ibin)) & (rad<binedges(ibin+1)));
        stampcb_ibin = stampcb(sp);
        stampps_ibin = stampps(sp);
        if spire
            [N,~] = histc(stampcb_ibin,Ibinedges_cb);
            N = N(1:end-1);  
        else
            [N,~] = histcounts(stampcb_ibin,Ibinedges_cb);
        end
        
        histcb(ibin,:) = histcb(ibin,:) + reshape(N,[1,numel(Ibinedges_cb)-1]);
        
        if spire
            [N,~] = histc(stampps_ibin,Ibinedges_ps);
            N = N(1:end-1);         
        else
            [N,~] = histcounts(stampps_ibin,Ibinedges_ps);
        end
        
        histps(ibin,:) = histps(ibin,:) + reshape(N,[1,numel(Ibinedges_ps)-1]);
        
    end
    
    if ismember(i,N_arr)
        Nidx = find(N_arr==i);
        histbkdat(Nidx).N = N_arr(Nidx);
        histbkdat(Nidx).histcb = histcb;
        histbkdat(Nidx).histps = histps;
    end
    
    
    if ismember(i,idx_print)
        print_count = print_count + 1;
        disp(sprintf('stack %d %% %d/%d random src'...
            ,print_count,i,max(N_arr)));         
    end
    
end

return