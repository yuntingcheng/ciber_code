function [stampercb,stamperps,hitmap]=...
    stackihl_sim0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
    mask_inst,strmask,strnum,unmask,Nsrc,verbose,isub,Nsub,idx_stack_arr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stack src based on sim catalog
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
dt=get_dark_times(flight,inst,ifield);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');
catfile=strcat(catdir,'sides.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
xg_arr=x_arr+1;
yg_arr=y_arr+1;

if inst==1
    mg_arr = squeeze(M(:,5)');
elseif inst==2
    mg_arr = squeeze(M(:,6)');
end

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/trilegal/');
catfile=strcat(catdir,dt.name, '.csv');
M = csvread(catfile,1);
x_arr=squeeze(M(:,3)');
y_arr=squeeze(M(:,4)');
xs_arr=x_arr+1;
ys_arr=y_arr+1;

if inst==1
    ms_arr = squeeze(M(:,1)');
elseif inst==2
    ms_arr = squeeze(M(:,2)');
end

x_arr = [xs_arr, xg_arr];
y_arr = [ys_arr, yg_arr];
m_arr = [ms_arr, mg_arr];
cls_arr = [-1*ones(size(ms_arr)), ones(size(mg_arr))];

%%% count the center pix map
xround_arr=round(x_arr);
yround_arr=round(y_arr);

centnum_map = zeros(720);
for i=1:numel(xround_arr)
    if xround_arr(i) < 720.5 & yround_arr(i)<720.5 & m_arr(i) < 20
    centnum_map(xround_arr(i),yround_arr(i))=...
        centnum_map(xround_arr(i),yround_arr(i))+1;
    end
end

%%% select sources %%%
sp=find(m_arr<=m_max & m_arr>m_min & cls_arr==type);
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

if numel(idx_stack_arr) == 0
    idx_stack_arr = idx_arr;
end

for i=idx_stack_arr
    cbmapi = cbmap.*strmask.*mask_inst;
    psmapi = psmap.*strmask.*mask_inst;
    %%% unmask the target
    if unmask
        radmap = make_radius_map(cbmapi,subx_arr(i),suby_arr(i));
        sp1 = find (radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
        cbmapi(sp1) = cbmap(sp1);
        psmapi(sp1) = psmap(sp1);
    end
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');    
    %%% rebin to 10x finer map
    stampcb0 = imresize(cbmapi,10,'nearest');
    stampps0 = imresize(psmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stampcb = stampcb0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stampps = stampps0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp = zeros(size(stampcb));
    maskstamp(find(stampcb)) = 1;

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
    if ismember(i,idx_print) & verbose>0
        print_count = print_count + 1;
        disp(sprintf('stack %d %% %d/%d src between %.1f<m<%.1f '...
            ,print_count*5,i,numel(idx_stack_arr),m_min,m_max));         
    end
end
%%

    disp(sprintf('stack %d src between %.1f<m<%.1f '...
        ,numel(idx_stack_arr),m_min,m_max));

return
