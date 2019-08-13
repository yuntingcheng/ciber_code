function [clipmax_arr, clipmin_arr,rbins]=...
stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
mask_inst,strnum,Nsrc,verbose,stackband,idx_stack_arr,rmin)
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
cls_arr(cls_arr==3)=1;
cls_arr(cls_arr==6)=-1;
photz_arr=squeeze(M(:,12)');

if strcmp(stackband,'y')
    m_arr=squeeze(M(:,9)');
elseif strcmp(stackband,'Ilin')
    m_arr=squeeze(M(:,14)');
elseif strcmp(stackband,'Hlin')
    m_arr=squeeze(M(:,15)');
elseif strcmp(stackband,'I')
    m_arr=squeeze(M(:,21)');
elseif strcmp(stackband,'H')
    m_arr=squeeze(M(:,22)');
end

sp=find(x_arr>0.5 & x_arr<1024.5 & y_arr>0.5 & y_arr<1024.5);

x_arr = x_arr(sp);
y_arr = y_arr(sp);
m_arr = m_arr(sp);
cls_arr = cls_arr(sp);
photz_arr=photz_arr(sp);

%%% count the center pix map
xround_arr=round(x_arr);
yround_arr=round(y_arr);

centnum_map = zeros(1024);
for i=1:numel(xround_arr)
    centnum_map(xround_arr(i),yround_arr(i))=...
        centnum_map(xround_arr(i),yround_arr(i))+1;
end

%%% select sources %%%
spg=find(m_arr<=m_max & m_arr>m_min & cls_arr==1 & photz_arr >= 0);
sps=find(m_arr<=m_max & m_arr>m_min & cls_arr==-1);
sp = [sps,spg];


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
