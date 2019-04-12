function [stamper,hitmap]=stackpsf_2m(flight,inst,ifield,m_min,m_max,dx,...
    cbmap,mask_inst,strmask,strnum,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stack src based on PanSTARRS catalog
%
%Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - type: 1:gal, -1:star
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - dx: stamp size is 2*dx+1, dx is 10 times finer CIBER pix
% - cbmap:
% - psmap:
% - mask_inst:
% - strmask:
% - strnum:
% - verbose:1 or 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PSC/');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

m_arr=squeeze(M(:,10)');% y band

sp=find(x_arr>0.5 & x_arr<1024.5 & y_arr>0.5 & y_arr<1024.5);

x_arr = x_arr(sp);
y_arr = y_arr(sp);
m_arr = m_arr(sp);

%%% count the center pix map
xround_arr=round(x_arr);
yround_arr=round(y_arr);

centnum_map = zeros(1024);
for i=1:numel(xround_arr)
    centnum_map(xround_arr(i),yround_arr(i))=...
        centnum_map(xround_arr(i),yround_arr(i))+1;
end

%%% select sources %%%
sp=find(m_arr<=m_max & m_arr>m_min);

submtot_arr=m_arr(sp);
subxtot_arr=x_arr(sp);
subytot_arr=y_arr(sp);

count_tot = numel(sp);

%%% source selection %%%
count_stack = 0;
subm_arr=[];
subx_arr=[];
suby_arr=[];
for i=1:numel(sp)
    use=true;
    
    xi = round(subxtot_arr(i));
    yi = round(subytot_arr(i));
    
    % center not coexist with other sources
    if centnum_map(xi,yi)~=1
        use = false;
    end
    
    % center not masked by sigclip mask
    if mask_inst(xi,yi)~=1
        use = false;
    end
    
    % center not masked by other sources
    if strnum(xi,yi)>1
        use = false;
    end
    
    % the center pixel is the max among 3x3 stamp
    cbmapi = cbmap.*strmask.*mask_inst;
    radi = get_mask_radius(inst,ifield,submtot_arr(i));% [arcsec]
    radmap = make_radius_map(cbmapi,subxtot_arr(i),subytot_arr(i));
    sp1 = find (radmap < radi./7 & strnum==1 & mask_inst==1);
    cbmapi(sp1) = cbmap(sp1);
    cbmapi = padarray(cbmapi,[1 1],0,'both');
    stamp = cbmapi(xi:xi+2,yi:yi+2);
    
    if find(stamp==max(stamp(:)))~=5
        use = false;
    end
    
    if use
        count_stack=count_stack+1;
        subm_arr=[subm_arr,submtot_arr(i)];
        subx_arr=[subx_arr,subxtot_arr(i)];
        suby_arr=[suby_arr,subytot_arr(i)];
    end
end

%%% set up stacking %%%
rad_arr = get_mask_radius(inst,ifield,subm_arr);% [arcsec]

stamper = zeros(2*dx+1);
hitmap = zeros(2*dx+1);

idx_arr = 1:numel(subm_arr);

%%% sig clip params %%%
nbins = 25;
iter_clip = 3;
sig = 3;
%%
for i=idx_arr
    cbmapi = cbmap.*strmask.*mask_inst;
    %%% unmask the target
    radmap = make_radius_map(cbmapi,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
    cbmapi(sp1) = cbmap(sp1);
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');    
    %%% rebin to 10x finer map
    stamp = imresize(cbmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stamp0 = stamp(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    
    %%% sigma clip
    rmin = 20; % don't clip with in rmin subpixels
    sig_clip_mask = stamp_clip(stamp0,dx+1,dx+1,nbins,sig,iter_clip,rmin);
    stamp = stamp0 .* sig_clip_mask;
        
    %%% stack
    stamper = stamper+stamp;
    hitmap = hitmap+sig_clip_mask;
    %%% print
    if verbose>0
        disp(sprintf('stack %d/%d(%d total) src between %.1f<m<%.1f '...
            ,i,count_stack,count_tot,m_min,m_max));    
    end

end

%%
    disp(sprintf('stack %d(%d total) src between %.1f<m<%.1f '...
        ,count_stack,count_tot,m_min,m_max));

return
