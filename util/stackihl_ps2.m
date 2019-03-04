function [stampercb,stamperps,hitmap,count_tot,count_stack]=...
    stackihl_ps2(flight,inst,ifield,type,m_min,m_max,dx,...
    cbmap,psmap,mask_inst,strmask,strnum,nbins,rmin,cliplim_arr,verbose)
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
%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PanSTARRS/');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

m_arr=squeeze(M(:,9)');
cls_arr=squeeze(M(:,10)');

sp=find(x_arr>0.5 & x_arr<1024.5 & y_arr>0.5 & y_arr<1024.5);

x_arr = x_arr(sp);
y_arr = y_arr(sp);
m_arr = m_arr(sp);
cls_arr = cls_arr(sp);

%%% count the center pix map
xround_arr=round(x_arr);
yround_arr=round(y_arr);

centnum_map = zeros(1024);
for i=1:numel(xround_arr)
    centnum_map(xround_arr(i),yround_arr(i))=...
        centnum_map(xround_arr(i),yround_arr(i))+1;
end

%%% select sources %%%
sp=find(m_arr<=m_max & m_arr>m_min & cls_arr==type);
if type ==0
    sp=find(m_arr<=m_max & m_arr>m_min);
end

if type==2
    sp=find(m_arr<=m_max & m_arr>m_min & cls_arr~=1 & cls_arr~=-1);
end

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
rad_arr = get_mask_radius(inst,ifield,subm_arr);

stampercb=zeros(2*dx+1);
stamperps=zeros(2*dx+1);
hitmap=zeros(2*dx+1);

idx_arr = 1:numel(subm_arr);

%%% remove the bad sources by hand %%%
if inst==1 & ifield==5 & m_min==18 & type==-1
    idx_arr(535)=[];
end
%%
if numel(idx_arr)>20
    idx_print=floor(numel(idx_arr)/20) * (1:20);
else
    idx_print=[];
end
print_count = 0;

for i=datasample(idx_arr,20,'Replace',false)%idx_arr
    cbmapi = cbmap.*strmask.*mask_inst;
    psmapi = psmap.*strmask.*mask_inst;
    %%% unmask the target
    radmap = make_radius_map(cbmapi,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
%     cbmapi(sp1) = cbmap(sp1);
%     psmapi(sp1) = psmap(sp1);
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');    
    %%% rebin to 10x finer map
    stampcb0 = imresize(cbmapi,10,'nearest');
    stampps0 = imresize(psmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
%     stampcb0 = stampcb0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
%     stampps0 = stampps0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
%     max_arr = cliplim_arr(1,:);
%     min_arr = cliplim_arr(2,:);
%     sig_clip_mask = stamp_clip2(stampcb0,max_arr,min_arr,dx+1,dx+1,nbins,rmin);
%     stampcb1 = stampcb0 .* sig_clip_mask;
%     stampps1 = stampps0 .* sig_clip_mask;
%     max_arr = cliplim_arr(3,:);
%     min_arr = cliplim_arr(4,:);
%     sig_clip_mask = stamp_clip2(stampps1,max_arr,min_arr,dx+1,dx+1,nbins,rmin);
%     stampcb = stampcb1 .* sig_clip_mask;
%     stampps = stampps1 .* sig_clip_mask;
%     maskstamp = sig_clip_mask;    
    stampcb = stampcb0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);%%
    stampps = stampps0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);%%
    maskstamp = zeros(size(stampcb));%%
    maskstamp(find(stampcb)) = 1;%%

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
    
    %%% print
%     if ismember(i,idx_print) & verbose>0
%         print_count = print_count + 1;
%         disp(sprintf('stack %d %% %d/%d(%d total) src between %.1f<m<%.1f '...
%             ,print_count*5,i,count_stack,count_tot,m_min,m_max));         
%     end
     print_count = print_count + 1;%%
    disp(sprintf('stack %d',print_count));%%
%     setwinsize(gcf,1000,400)
%     subplot(121)
%     imageclip(stampercb(dx-100:dx+100,dx-100:dx+100));
%     subplot(122)
%     imageclip(stamperps(dx-100:dx+100,dx-100:dx+100));
%     title(i)
%     drawnow
%     pause
    
end
%%

    disp(sprintf('stack %d(%d total) src between %.1f<m<%.1f '...
        ,count_stack,count_tot,m_min,m_max));

return
