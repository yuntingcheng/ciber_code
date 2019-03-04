function [stampercb,stamperps,maskstamper,count_tot,count_stack]=...
    stackihl_sides(flight,inst,ifield,m_min,m_max,dx,...
    cbmap,psmap,mask_inst,strmask,strnum,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stack src based on PanSTARRS catalog
%
%Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - dx: stamp size is 2*dx+1, dx is 10 times finer CIBER pix
% - cbmap:
% - psmap:
% - mask_inst:
% - strmask:
% - strnum:
% - verbose:1 or 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%

catfile=strcat(catdir,'sides.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

m_arr=squeeze(M(:,4)');

sp=find(m_arr<20);

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

%%% select sources not coexist with others in the same pixel %%%
count_stack = 0;
subm_arr=[];
subx_arr=[];
suby_arr=[];
for i=1:numel(sp)
    if centnum_map(round(subxtot_arr(i)),round(subytot_arr(i)))==1 ...        
            & mask_inst(round(subxtot_arr(i)),round(subytot_arr(i)))==1
        count_stack=count_stack+1;
        subm_arr=[subm_arr,submtot_arr(i)];
        subx_arr=[subx_arr,subxtot_arr(i)];
        suby_arr=[suby_arr,subytot_arr(i)];
    end
end

%%% set up stacking %%%
rad_arr = -3.08e4.*exp(-(subm_arr-10.84).^2./(5.253).^2) + ...
    3.091e4.*exp(-(subm_arr-10.83).^2./(5.26).^2); % [arcsec]
rad_arr(find(rad_arr<7))=7;

stampercb=zeros(2*dx+1);
stamperps=zeros(2*dx+1);
maskstamper=zeros(2*dx+1);

idx_arr = 1:numel(subm_arr);

%%% sig clip params %%%
nbins = 25;
iter_clip = 3;
sig = 3;

for i=idx_arr
    cbmapi = cbmap.*strmask.*mask_inst;
    psmapi = psmap.*strmask.*mask_inst;
    %%% unmask the target
    radmap = make_radius_map(cbmapi,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
    cbmapi(sp1) = cbmap(sp1);
    psmapi(sp1) = psmap(sp1);
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');
    %%% sigma clip
    xcent = round(subx_arr(i)) + dx/10;
    ycent = round(suby_arr(i)) + dx/10;
    mask_clip = stamp_clip(cbmapi,xcent,ycent,nbins,sig,iter_clip);
    %disp(sprintf('%d,%d',i,numel(find(cbmapi))-numel(find(mask_clip))));
    cbmapi = cbmapi.*mask_clip;
    psmapi = psmapi.*mask_clip;
    %%% rebin to 10x finer map
    stampcb = imresize(cbmapi,10,'nearest');
    stampps = imresize(psmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stampcb = stampcb(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stampps = stampps(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp = zeros(size(stampcb));
    maskstamp(find(stampcb)) = 1;
    
    %%% stack
    stampercb=stampercb+stampcb;
    stamperps=stamperps+stampps;
    maskstamper=maskstamper+maskstamp;
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
