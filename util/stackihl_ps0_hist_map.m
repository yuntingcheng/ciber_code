function [stampercb,stamperps,hitmap,...
    profcb_arr,profps_arr,hitmap_arr]=...
stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,strmask,...
strnum,unmask,verbose,rmin,clipmaxs,clipmins,x_arr,y_arr,m_arr,subpix)
%%

if isnan(rmin)
    rad_arr = get_mask_radius_th(inst,ifield,m_arr,0.5);
else
    rad_arr = get_mask_radius_th(inst,ifield,m_arr,0.5,'rmin',rmin);
end


idx_arr = 1:numel(m_arr);

if numel(idx_arr)>20
    idx_print=floor(numel(idx_arr)/20) * (1:20);
else
    idx_print=[];
end
print_count = 0;

nbins = 25;
%%
if subpix
    stampercb=zeros(2*dx+1);
    stamperps=zeros(2*dx+1);
    hitmap=zeros(2*dx+1);
    rad = make_radius_map(zeros(2*dx+1),dx+1,dx+1);
    profile = radial_prof(rad,ones(2*dx+1),dx+1,dx+1,1,nbins);
    rbinedges = profile.binedges;
else
    stampercb=zeros(2*dx/10+1);
    stamperps=zeros(2*dx/10+1);
    hitmap=zeros(2*dx/10+1);
    rad = make_radius_map(zeros(2*dx+1),dx+1,dx+1);
    profile = radial_prof(rad,ones(2*dx+1),dx+1,dx+1,1,nbins);
    rbinedges = profile.binedges./10;
    rad = make_radius_map(zeros(2*dx/10+1),dx/10+1,dx/10+1);
end

profcb_arr = zeros([1,nbins]);
profps_arr = zeros([1,nbins]);
hitmap_arr = zeros([1,nbins]);

% fig = figure;
% setwinsize(gcf,800,300);
% gifname=strcat('/Users/ytcheng/Downloads/stack.gif');
for i=idx_arr
    cbmapi = cbmap.*strmask.*mask_inst;
    psmapi = psmap.*strmask.*mask_inst;
    
    %%% unmask the target and clip
    if unmask
        radmap = make_radius_map(cbmapi,x_arr(i),y_arr(i));
        sp1 = find (radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
        cbmapi(sp1) = cbmap(sp1);
        psmapi(sp1) = psmap(sp1);
        unmaskpix = zeros(size(cbmap));
        unmaskpix(sp1)=1;
        if numel(sp1)~=0
            for ibin=1:nbins
                spi = find (unmaskpix==1 & radmap.*10>=rbinedges(ibin) ...
                    & radmap.*10<rbinedges(ibin+1) & cbmap>clipmaxs(1,ibin));
                cbmapi(spi)=0;psmapi(spi)=0;
                spi = find (unmaskpix==1 & radmap.*10>=rbinedges(ibin) ...
                    & radmap.*10<rbinedges(ibin+1) & cbmap<clipmins(1,ibin));
                cbmapi(spi)=0;psmapi(spi)=0;
                spi = find (unmaskpix==1 & radmap.*10>=rbinedges(ibin) ...
                    & radmap.*10<rbinedges(ibin+1) & psmap>clipmaxs(2,ibin));
                cbmapi(spi)=0;psmapi(spi)=0;
                spi = find (unmaskpix==1 & radmap.*10>=rbinedges(ibin) ...
                    & radmap.*10<rbinedges(ibin+1) & psmap<clipmins(2,ibin));
                cbmapi(spi)=0;psmapi(spi)=0;                
            end 
        end
    end
    
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
        sp = find((maskstamp~=0) & (rad>=rbinedges(ibin))...
            & (rad<rbinedges(ibin+1)));
        profcb_arr(ibin) = profcb_arr(ibin) + sum(stampcb(sp));
        profps_arr(ibin) = profps_arr(ibin) + sum(stampps(sp));
        hitmap_arr(ibin) = hitmap_arr(ibin) + numel(sp);
    end

    stampcb = stampcb.*maskstamp;
    stampps = stampps.*maskstamp;
    %%% stack
    if mod(i,4)==0
        stampercb=stampercb+imrotate(stampcb, 0);
        stamperps=stamperps+imrotate(stampps, 0);
        hitmap=hitmap+imrotate(maskstamp, 0);
    elseif mod(i,4)==1
        stampercb=stampercb+imrotate(stampcb, 90);
        stamperps=stamperps+imrotate(stampps, 90);
        hitmap=hitmap+imrotate(maskstamp, 90);
    elseif mod(i,4)==2
        stampercb=stampercb+imrotate(stampcb, 180);
        stamperps=stamperps+imrotate(stampps, 180);
        hitmap=hitmap+imrotate(maskstamp, 180);
    elseif mod(i,4)==3
        stampercb=stampercb+imrotate(stampcb, 270);
        stamperps=stamperps+imrotate(stampps, 270);
        hitmap=hitmap+imrotate(maskstamp, 270);
    end

    
%     % make gif
%     if numel(find(maskstamp)~=0)>0
%         subplot(1,2,1)
%         imageclip(stampercb./hitmap);
%         colorbar('off');
%         subplot(1,2,2)
%         imageclip(stamperps./hitmap);
%         colorbar('off');
%         frame = getframe(fig);
%         [A,map] = rgb2ind(frame2im(frame),256);
%         if i == 1
%             imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.2);
%         else
%             imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.2);
%         end
%     end
    %%% print
    if ismember(i,idx_print) & verbose
        print_count = print_count + 1;
        disp(sprintf('stack %d %% %d/%d src',print_count*5,i,numel(idx_arr)));
    end
end
%%
profcb_arr = profcb_arr./hitmap_arr;
profps_arr = profps_arr./hitmap_arr;

% disp(sprintf('stack %d src',numel(idx_arr)));
return
