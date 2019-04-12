function [stampercb,stamperps,hitmap,idx_stack_arr]=...
    stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,dx,...
    cbmap,psmap,mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
    cliplim_arrout,binedgesout,Nsrc,verbose,isub,Nsub,stackband)
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

cls_arr=squeeze(M(:,11)');
if type == 1
    type = 3;
elseif type == -1
    type =6;
end
    

if strcmp(stackband,'y')
    m_arr=squeeze(M(:,9)');
elseif strcmp(stackband,'Ilin')
    m_arr=squeeze(M(:,14)');
elseif strcmp(stackband,'Hlin')
    m_arr=squeeze(M(:,15)');
elseif strcmp(stackband,'I')
    m_arr=squeeze(M(:,16)');
elseif strcmp(stackband,'H')
    m_arr=squeeze(M(:,17)');
end


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
%     sp=find(m_arr<=m_max & m_arr>m_min & cls_arr~=1 & cls_arr~=-1);
    sp=find(m_arr<=m_max & m_arr>m_min & cls_arr~=3 & cls_arr~=6);
end

submtot_arr=m_arr(sp);
subxtot_arr=x_arr(sp);
subytot_arr=y_arr(sp);

%%% select sources not coexist with others in the same pixel %%%
subm_arr=[];
subx_arr=[];
suby_arr=[];
for i=1:numel(sp)
    if centnum_map(round(subxtot_arr(i)),round(subytot_arr(i)))==1 ...        
            & mask_inst(round(subxtot_arr(i)),round(subytot_arr(i)))==1
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


if numel(idx_arr)>20
    idx_print=floor(numel(idx_arr)/20) * (1:20);
else
    idx_print=[];
end
print_count = 0;

if Nsrc ~= 0 & Nsrc < numel(idx_arr)
    idx_arr = datasample(idx_arr,Nsrc,'Replace',false);
end

if Nsub ~= 0
    idx_arr = idx_arr(isub:Nsub:floor(numel(idx_arr)/Nsub)*Nsub);
end
idx_stack_arr = [];
for i=idx_arr
    cbmapi = cbmap.*strmask.*mask_inst;
    psmapi = psmap.*strmask.*mask_inst;
    %%% unmask the target
    radmap = make_radius_map(cbmapi,subx_arr(i),suby_arr(i));
    sp1 = find(radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
    cbmapi(sp1) = cbmap(sp1);
    psmapi(sp1) = psmap(sp1);
    targmask = zeros(1024);
    targmask(sp1) = 1;
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');
    targmask = padarray(targmask,[dx/10 dx/10],0,'both');
    %%% rebin to 10x finer map
    stampcb0 = imresize(cbmapi,10,'nearest');
    stampps0 = imresize(psmapi,10,'nearest');
    targmask = imresize(targmask,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stampcb = stampcb0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stampps = stampps0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    targmask = targmask(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp = zeros(size(stampcb));
    maskstamp(find(stampcb)) = 1;

    %%% clip
    rad = make_radius_map(stampcb,dx+1,dx+1);
    for ibin=1:numel(binedgesout) - 1
        clipmaxcb = cliplim_arrout(1,ibin);
        clipmincb = cliplim_arrout(2,ibin);
        clipmaxps = cliplim_arrout(3,ibin);
        clipminps = cliplim_arrout(4,ibin);
        
        spclip = find(rad >= binedgesout(ibin) ...
        & rad < binedgesout(ibin+1) & maskstamp & targmask ...
        & (stampcb >= clipmaxcb | stampcb <= clipmincb |...
        stampps >= clipmaxps | stampps <= clipminps));
        clip_val = unique(stampcb(spclip))';
        maskstamp(ismember(stampcb,clip_val)) = 0;
    end
    
    %%% clip the whole source if systematically brighter all bins
    bin_sig_cb = [];
    bin_sig_ps = [];
    for ibin=1:numel(binedgesin) - 1
        clipmucb = (cliplim_arrin(1,ibin) + cliplim_arrin(2,ibin))./2;
        clipsigcb = (cliplim_arrin(1,ibin) - cliplim_arrin(2,ibin))./2;
        clipmups = (cliplim_arrin(3,ibin) + cliplim_arrin(4,ibin))./2;
        clipsigps = (cliplim_arrin(3,ibin) - cliplim_arrin(4,ibin))./2;

        spclip = find(rad >= binedgesin(ibin) ...
        & rad < binedgesin(ibin+1) & maskstamp);
    
        bin_sig_cb = [bin_sig_cb,(mean(stampcb(spclip)) - clipmucb)./clipsigcb];
        bin_sig_ps = [bin_sig_ps,(mean(stampps(spclip)) - clipmups)./clipsigps];
        
    end
    
    bin_sig_cb = bin_sig_cb(1:end-1);
    bin_sig_ps = bin_sig_ps(1:end-1);
    if all(bin_sig_cb > 20) | all(bin_sig_ps > 5)
        maskstamp = maskstamp.*0;
    else
        idx_stack_arr = [idx_stack_arr,i];

    end
    
    stampps = stampps.*maskstamp;
    stampcb = stampcb.*maskstamp;
    
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
    
%     figure
%     imageclip(stampercb./hitmap);
%     title(i);
%     drawnow

    %%% print
    if ismember(i,idx_print) & verbose>0
        print_count = print_count + 1;
        disp(sprintf('stack %d %% %d/%d src between %.1f<m<%.1f '...
            ,print_count*5,i,numel(idx_arr),m_min,m_max));         
    end
    
end
%%

    disp(sprintf('stack %d src between %.1f<m<%.1f '...
        ,numel(idx_stack_arr),m_min,m_max));

return
