function [srcdat]=ps_src_select(flight,inst,ifield,m_min,m_max,...
    mask_inst,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select sources to be stack
% sample_type:
% - 'all': all the sources
% - 'jack_random': randomly separate sources into Nsub subsets
% - 'jack_region': separate sources into 16(4x4) spatial regions
% - 'jack_random_uni': randomly separate uniformizec sources into Nsub subsets
% - 'jack_region_uni': separate uniformized sources into 16(4x4) spatial regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('m_min',@isnumeric);
  p.addRequired('m_max',@isnumeric);
  p.addRequired('mask_inst',@isnumeric);
  p.addOptional('Nsub',16,@isnumeric);
  p.addOptional('sample_type','all',@ischar);
  p.addOptional('make_plot',false,@islogical);
  p.addOptional('HSC',false,@islogical);
  p.parse(flight,inst,ifield,m_min,m_max,mask_inst,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  m_min    = p.Results.m_min;
  m_max    = p.Results.m_max;
  mask_inst= p.Results.mask_inst;
  sample_type = p.Results.sample_type;
  Nsub     = p.Results.Nsub;
  make_plot= p.Results.make_plot;
  HSC      = p.Results.HSC;
  clear p varargin;
%%
mypaths=get_paths(flight);

if ~HSC
    catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/PanSTARRS/');

    dt=get_dark_times(flight,inst,ifield);

    %%% read cat data %%%
    catfile=strcat(catdir,dt.name,'.csv');

    M = csvread(catfile,1);
    x1_arr=squeeze(M(:,4)');
    y1_arr=squeeze(M(:,3)');
    x2_arr=squeeze(M(:,6)');
    y2_arr=squeeze(M(:,5)');    

    x1_arr=x1_arr+1;
    y1_arr=y1_arr+1;
    x2_arr=x2_arr+1;
    y2_arr=y2_arr+1;

    cls_arr=squeeze(M(:,13)');
    cls_arr(cls_arr==3)=1;
    cls_arr(cls_arr==6)=-1;
    photz_arr=squeeze(M(:,14)');

%     if strcmp(stackband,'y')
%         m_arr=squeeze(M(:,11)');
%     elseif strcmp(stackband,'Ilin')
%         m_arr=squeeze(M(:,16)');
%     elseif strcmp(stackband,'Hlin')
%         m_arr=squeeze(M(:,17)');
%     elseif strcmp(stackband,'I')
%         m_arr=squeeze(M(:,23)');
%     elseif strcmp(stackband,'H')
%         m_arr=squeeze(M(:,24)');
%     end
    m_arr=squeeze(M(:,23)');% always select by I band
    
    sp=find(x1_arr>0.5 & x1_arr<1024.5 & y1_arr>0.5 & y1_arr<1024.5 & ...
        x2_arr>0.5 & x2_arr<1024.5 & y2_arr>0.5 & y2_arr<1024.5);

    x1_arr = x1_arr(sp);
    y1_arr = y1_arr(sp);
    x2_arr = x2_arr(sp);
    y2_arr = y2_arr(sp);
    m_arr = m_arr(sp);
    cls_arr = cls_arr(sp);
    photz_arr=photz_arr(sp);

    %%% count the center pix map
    x1round_arr=round(x1_arr);
    y1round_arr=round(y1_arr);
    x2round_arr=round(x2_arr);
    y2round_arr=round(y2_arr);
    centnum_map1 = zeros(1024);
    centnum_map2 = zeros(1024);
    for i=1:numel(x1round_arr)
        centnum_map1(x1round_arr(i),y1round_arr(i))=...
            centnum_map1(x1round_arr(i),y1round_arr(i))+1;
    end
    for i=1:numel(x2round_arr)
        centnum_map2(x2round_arr(i),y2round_arr(i))=...
            centnum_map2(x2round_arr(i),y2round_arr(i))+1;
    end
    
    spg=find(m_arr<=m_max & m_arr>m_min & cls_arr==1 & photz_arr >= 0);
    sps=find(m_arr<=m_max & m_arr>m_min & cls_arr==-1);
    sp = [sps,spg];

    m_arr=m_arr(sp);
    x1_arr=x1_arr(sp);
    y1_arr=y1_arr(sp);
    x2_arr=x2_arr(sp);
    y2_arr=y2_arr(sp);
    z_arr=photz_arr(sp);
    cls_arr=ones(size(m_arr));
    cls_arr(1:numel(sps))=-1;

    %%% select sources not coexist with others in the same pixel %%%
    mask_inst1 = squeeze(mask_inst(1,:,:));
    mask_inst2 = squeeze(mask_inst(2,:,:));
    subm_arr=[];
    subx1_arr=[];
    suby1_arr=[];
    subx2_arr=[];
    suby2_arr=[];
    subz_arr=[];
    subcls_arr=[];
    for i=1:numel(sp)
        if centnum_map1(round(x1_arr(i)),round(y1_arr(i)))==1 ...
                & centnum_map2(round(x2_arr(i)),round(y2_arr(i)))==1 ...
                & mask_inst1(round(x1_arr(i)),round(y1_arr(i)))==1 ...
                & mask_inst2(round(x2_arr(i)),round(y2_arr(i)))==1
            subm_arr=[subm_arr,m_arr(i)];
            subx1_arr=[subx1_arr,x1_arr(i)];
            suby1_arr=[suby1_arr,y1_arr(i)];
            subx2_arr=[subx2_arr,x2_arr(i)];
            suby2_arr=[suby2_arr,y2_arr(i)];
            subz_arr=[subz_arr,z_arr(i)];
            subcls_arr=[subcls_arr,cls_arr(i)];
        end
    end
    
    randidx = randperm(numel(subm_arr));
    if inst==1
        x_arr = subx1_arr(randidx);
        y_arr = suby1_arr(randidx);
    else
        x_arr = subx2_arr(randidx);
        y_arr = suby2_arr(randidx);        
    end
    
    z_arr = subz_arr(randidx);
    m_arr = subm_arr(randidx);
    cls_arr = subcls_arr(randidx);

    xg_arr = x_arr(cls_arr==1);
    yg_arr = y_arr(cls_arr==1);
    mg_arr = m_arr(cls_arr==1);
    zg_arr = z_arr(cls_arr==1);
    xs_arr = x_arr(cls_arr==-1);
    ys_arr = y_arr(cls_arr==-1);
    ms_arr = m_arr(cls_arr==-1);
   
else
    hsc_idx = ifield;
    catdir=mypaths.hsccatdir;
    [field,mask_inst] = HSC_fields_info(hsc_idx);

    catfile=strcat(catdir,field,'.csv');

    M = csvread(catfile,1);

    x_arr=squeeze(M(:,3)');
    y_arr=squeeze(M(:,4)');
    x_arr=x_arr-1.5;
    y_arr=y_arr-1.5;
    
    cls_arr=squeeze(M(:,12)');
    photz_arr=squeeze(M(:,14)');

%     if inst==1
%         m_arr=squeeze(M(:,10)');
%     else
%         m_arr=squeeze(M(:,11)');
%     end
    m_arr=squeeze(M(:,10)');% always select by I band mag
    
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
        if m_arr(i)<20
            centnum_map(xround_arr(i),yround_arr(i))=...
                centnum_map(xround_arr(i),yround_arr(i))+1;
        end
    end
    spg=find(m_arr<=m_max & m_arr>m_min & cls_arr==1 & photz_arr >= 0);
    sps=find(m_arr<=m_max & m_arr>m_min & cls_arr==-1);
    sp = [sps,spg];

    m_arr=m_arr(sp);
    x_arr=x_arr(sp);
    y_arr=y_arr(sp);
    z_arr=photz_arr(sp);
    cls_arr=ones(size(x_arr));
    cls_arr(1:numel(sps))=-1;

    %%% select sources not coexist with others in the same pixel %%%
    subm_arr=[];
    subx_arr=[];
    suby_arr=[];
    subz_arr=[];
    subcls_arr=[];
    for i=1:numel(sp)
        if centnum_map(round(x_arr(i)),round(y_arr(i)))==1 ...        
                & mask_inst(round(x_arr(i)),round(y_arr(i)))==1
            subm_arr=[subm_arr,m_arr(i)];
            subx_arr=[subx_arr,x_arr(i)];
            suby_arr=[suby_arr,y_arr(i)];
            subz_arr=[subz_arr,z_arr(i)];
            subcls_arr=[subcls_arr,cls_arr(i)];
        end
    end
    
    randidx = randperm(numel(subm_arr));
    x_arr = subx_arr(randidx);
    y_arr = suby_arr(randidx);
    z_arr = subz_arr(randidx);
    m_arr = subm_arr(randidx);
    cls_arr = subcls_arr(randidx);

    xg_arr = x_arr(cls_arr==1);
    yg_arr = y_arr(cls_arr==1);
    mg_arr = m_arr(cls_arr==1);
    zg_arr = z_arr(cls_arr==1);
    xs_arr = x_arr(cls_arr==-1);
    ys_arr = y_arr(cls_arr==-1);
    ms_arr = m_arr(cls_arr==-1);
end
%%
if strcmp(sample_type,'all')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ng = numel(xg_arr);
    srcdat.Ns = numel(xs_arr);
    srcdat.xg_arr = xg_arr;
    srcdat.yg_arr = yg_arr;
    srcdat.mg_arr = mg_arr;
    srcdat.zg_arr = zg_arr;
    srcdat.xs_arr = xs_arr;
    srcdat.ys_arr = ys_arr;
    srcdat.ms_arr = ms_arr;    
    return
end
%%
if strcmp(sample_type,'jack_random')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ng = numel(xg_arr);
    srcdat.Ns = numel(xs_arr);
    for isub=1:Nsub
        srcdat.sub(isub).xg_arr = xg_arr(isub:Nsub:numel(xg_arr));
        srcdat.sub(isub).yg_arr = yg_arr(isub:Nsub:numel(xg_arr));
        srcdat.sub(isub).mg_arr = mg_arr(isub:Nsub:numel(xg_arr));
        srcdat.sub(isub).zg_arr = zg_arr(isub:Nsub:numel(xg_arr));
        srcdat.sub(isub).Ng = numel(srcdat.sub(isub).xg_arr);
        srcdat.sub(isub).xs_arr = xs_arr(isub:Nsub:numel(xs_arr));
        srcdat.sub(isub).ys_arr = ys_arr(isub:Nsub:numel(xs_arr));
        srcdat.sub(isub).ms_arr = ms_arr(isub:Nsub:numel(xs_arr));
        srcdat.sub(isub).Ns = numel(srcdat.sub(isub).xs_arr);
    end
    return
end
%%
if strcmp(sample_type,'jack_region')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ng = numel(xg_arr);
    srcdat.Ns = numel(xs_arr);
    isub=0;
    for i=1:4
        for j=1:4
            isub = isub + 1;
            spg = find((xg_arr>=(i*256-255.5)) & (xg_arr<(i*256+0.5)) ...
                & (yg_arr>=(j*256-255.5)) & (yg_arr<(j*256.5)));
            srcdat.sub(isub).xg_arr = xg_arr(spg);
            srcdat.sub(isub).yg_arr = yg_arr(spg);
            srcdat.sub(isub).mg_arr = mg_arr(spg);
            srcdat.sub(isub).zg_arr = zg_arr(spg);
            srcdat.sub(isub).Ng = numel(srcdat.sub(isub).xg_arr);
            sps = find((xs_arr>=(i*256-255.5)) & (xs_arr<(i*256+0.5)) ...
                & (ys_arr>=(j*256-255.5)) & (ys_arr<(j*256.5)));
            srcdat.sub(isub).xs_arr = xs_arr(sps);
            srcdat.sub(isub).ys_arr = ys_arr(sps);
            srcdat.sub(isub).ms_arr = ms_arr(sps);
            srcdat.sub(isub).Ns = numel(srcdat.sub(isub).xs_arr);               
        end
    end
end
%%
%%% select uniform sources %%%
hitmaps = zeros(1024);
hitmapg = zeros(1024);
for i=1:numel(xs_arr)
    hitmaps(round(xs_arr(i)),round(ys_arr(i))) = ...
        hitmaps(round(xs_arr(i)),round(ys_arr(i))) + 1;
end
hitmaps = hitmaps.*mask_inst;
for i=1:numel(xg_arr)
    hitmapg(round(xg_arr(i)),round(yg_arr(i))) = ...
        hitmapg(round(xg_arr(i)),round(yg_arr(i))) + 1;
end
hitmapg = hitmapg.*mask_inst;

Nreg = 1;
binmapg = rebin_map_coarse(hitmapg,1024/Nreg);
binmaps = rebin_map_coarse(hitmaps,1024/Nreg);
binmapg(binmapg>0)=1;
binmaps(binmaps>0)=1;
while numel(find(binmapg~=binmaps))/numel(binmapg) < 0.1
    Nreg = Nreg*2;
    binmapg = rebin_map_coarse(hitmapg,1024/Nreg).*((1024/Nreg).^2);
    binmaps = rebin_map_coarse(hitmaps,1024/Nreg).*((1024/Nreg).^2);
    binmapg(binmapg>0)=1;
    binmaps(binmaps>0)=1;
end
Nreg = Nreg/2;
Npixsub = 1024/Nreg;
binmapg = rebin_map_coarse(hitmapg,1024/Nreg).*((1024/Nreg).^2);
binmaps = rebin_map_coarse(hitmaps,1024/Nreg).*((1024/Nreg).^2);
binmapg(binmapg>0)=1;
binmaps(binmaps>0)=1;

binmask = binmapg.*binmaps;
binmapg = rebin_map_coarse(hitmapg,1024/Nreg).*((1024/Nreg).^2);
binmaps = rebin_map_coarse(hitmaps,1024/Nreg).*((1024/Nreg).^2);
Nsrc_subs = min(binmaps(binmask~=0));
Nsrc_subg = min(binmapg(binmask~=0));

subxg_arr = [];
subyg_arr = [];
submg_arr = [];
subzg_arr = [];
subxs_arr = [];
subys_arr = [];
subms_arr = [];

for i=1:Nreg
    for j=1:Nreg
        if binmask(i,j)==0
            continue
        end
        spg = find((round(xg_arr)>=(i*Npixsub-(Npixsub-1))) ...
            & (round(xg_arr)<=(i*Npixsub)) ...
            & (round(yg_arr)>=(j*Npixsub-(Npixsub-1))) ...
            & (round(yg_arr)<=(j*Npixsub)));
        sps = find((round(xs_arr)>=(i*Npixsub-(Npixsub-1))) ...
            & (round(xs_arr)<=(i*Npixsub)) ...
            & (round(ys_arr)>=(j*Npixsub-(Npixsub-1))) ...
            & (round(ys_arr)<=(j*Npixsub)));
        
        if Nsrc_subs~=1 | Nsrc_subs~=1
            spg = spg(randsample(numel(spg),Nsrc_subg));  
            sps = sps(randsample(numel(sps),Nsrc_subs));
        else
            Nsrcij = min(numel(spg),numel(sps));
            spg = spg(randsample(numel(spg),Nsrcij));
            sps = sps(randsample(numel(sps),Nsrcij));
        end
        
        subxg_arr = [subxg_arr, xg_arr(spg)];
        subyg_arr = [subyg_arr, yg_arr(spg)];
        subzg_arr = [subzg_arr, zg_arr(spg)];
        submg_arr = [submg_arr, mg_arr(spg)];

        
        subxs_arr = [subxs_arr, xs_arr(sps)];
        subys_arr = [subys_arr, ys_arr(sps)];
        subms_arr = [subms_arr, ms_arr(sps)];

    end
end
dropfrac = 0.75;
randidx = randperm(numel(subxs_arr));
randidx = randidx(1:round(numel(randidx)*dropfrac));
subxs_arr = subxs_arr(randidx);
subys_arr = subys_arr(randidx);
subms_arr = subms_arr(randidx);
randidx = randperm(numel(subxg_arr));
randidx = randidx(1:round(numel(randidx)*dropfrac));
subxg_arr = subxg_arr(randidx);
subyg_arr = subyg_arr(randidx);
submg_arr = submg_arr(randidx);
subzg_arr = subzg_arr(randidx);
%% plot the source before and after selection
if make_plot
    figure
    setwinsize(gcf,1000,1000)

    hitmaps = zeros(1024);
    hitmapg = zeros(1024);
    for i=1:numel(xs_arr)
        hitmaps(round(xs_arr(i)),round(ys_arr(i))) = ...
            hitmaps(round(xs_arr(i)),round(ys_arr(i))) + 1;
    end
    hitmaps = hitmaps.*mask_inst;
    for i=1:numel(xg_arr)
        hitmapg(round(xg_arr(i)),round(yg_arr(i))) = ...
            hitmapg(round(xg_arr(i)),round(yg_arr(i))) + 1;
    end
    hitmapg = hitmapg.*mask_inst;
    binmapg = rebin_map_coarse(hitmapg,8).*8^2;
    binmaps = rebin_map_coarse(hitmaps,8).*8^2;
    subplot(2,2,1)
    imageclip(binmaps);
    title(strcat(num2str(numel(xs_arr)),' stars'),'fontsize',15);

    caxis([0,3]);
    subplot(2,2,2)
    imageclip(binmapg);
    title(strcat(num2str(numel(xg_arr)),' galaxies'),'fontsize',15);
    caxis([0,3]);

    hitmaps = zeros(1024);
    hitmapg = zeros(1024);
    for i=1:numel(subxs_arr)
        hitmaps(round(subxs_arr(i)),round(subys_arr(i))) = ...
            hitmaps(round(subxs_arr(i)),round(subys_arr(i))) + 1;
    end
    hitmaps = hitmaps.*mask_inst;
    for i=1:numel(subxg_arr)
        hitmapg(round(subxg_arr(i)),round(subyg_arr(i))) = ...
            hitmapg(round(subxg_arr(i)),round(subyg_arr(i))) + 1;
    end
    hitmapg = hitmapg.*mask_inst;
    binmapg = rebin_map_coarse(hitmapg,8).*8^2;
    binmaps = rebin_map_coarse(hitmaps,8).*8^2;
    subplot(2,2,3)
    imageclip(binmaps);
    title(strcat(num2str(numel(subxs_arr)),' stars'),'fontsize',15);
    caxis([0,3]);
    subplot(2,2,4)
    imageclip(binmapg);
    title(strcat(num2str(numel(subxg_arr)),' galaxies'),'fontsize',15);
    caxis([0,3]);
    s=suptitle(strcat(dt.name,{'   '},num2str(m_min),'<m<',num2str(m_max)));
    set(s,'FontSize',30);
end
%%
if strcmp(sample_type,'jack_random_uni')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ngtot = numel(xg_arr);
    srcdat.Nstot = numel(xs_arr);
    srcdat.Ng = numel(subxg_arr);
    srcdat.Ns = numel(subxs_arr);
    srcdat.Nreg_uni = Nreg;
    srcdat.dropfrac = dropfrac;
    for isub=1:Nsub
        srcdat.sub(isub).xg_arr = subxg_arr(isub:Nsub:numel(subxg_arr));
        srcdat.sub(isub).yg_arr = subyg_arr(isub:Nsub:numel(subxg_arr));
        srcdat.sub(isub).mg_arr = submg_arr(isub:Nsub:numel(subxg_arr));
        srcdat.sub(isub).zg_arr = subzg_arr(isub:Nsub:numel(subxg_arr));
        srcdat.sub(isub).Ng = numel(srcdat.sub(isub).xg_arr);
        srcdat.sub(isub).xs_arr = subxs_arr(isub:Nsub:numel(subxs_arr));
        srcdat.sub(isub).ys_arr = subys_arr(isub:Nsub:numel(subxs_arr));
        srcdat.sub(isub).ms_arr = subms_arr(isub:Nsub:numel(subxs_arr));
        srcdat.sub(isub).Ns = numel(srcdat.sub(isub).xs_arr);
    end
    return
end
%%
if strcmp(sample_type,'jack_region_uni')
    srcdat.sample_type = sample_type;
    srcdat.m_min = m_min;
    srcdat.m_max = m_max;
    srcdat.Ngtot = numel(xg_arr);
    srcdat.Nstot = numel(xs_arr);
    srcdat.Ng = numel(subxg_arr);
    srcdat.Ns = numel(subxs_arr);
    srcdat.Nreg_uni = Nreg;
    srcdat.dropfrac = dropfrac;
    isub=0;
    for i=1:4
        for j=1:4
            isub = isub + 1;
            spg = find((subxg_arr>=(i*256-255.5)) & (subxg_arr<(i*256+0.5)) ...
                & (subyg_arr>=(j*256-255.5)) & (subyg_arr<(j*256.5)));
            srcdat.sub(isub).xg_arr = subxg_arr(spg);
            srcdat.sub(isub).yg_arr = subyg_arr(spg);
            srcdat.sub(isub).mg_arr = submg_arr(spg);
            srcdat.sub(isub).zg_arr = subzg_arr(spg);
            srcdat.sub(isub).Ng = numel(srcdat.sub(isub).xg_arr);
            sps = find((subxs_arr>=(i*256-255.5)) & (subxs_arr<(i*256+0.5)) ...
                & (subys_arr>=(j*256-255.5)) & (subys_arr<(j*256.5)));
            srcdat.sub(isub).xs_arr = subxs_arr(sps);
            srcdat.sub(isub).ys_arr = subys_arr(sps);
            srcdat.sub(isub).ms_arr = subms_arr(sps);
            srcdat.sub(isub).Ns = numel(srcdat.sub(isub).xs_arr);               
        end
    end
    return
end
%%

return
