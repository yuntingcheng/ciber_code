flight = 40030;
inst = 2;

for ifield = 4:8

mypaths=get_paths(flight);

m_min_arr = [0,8:22];
m_max_arr = [8:23];

m_min_arr_2m = [0,7:21];
m_max_arr_2m = [7:22];

dt=get_dark_times(flight,inst,ifield);

%%% get srcmap %%%
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

map = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps1_all.fits'));

%%% get cats %%%
catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PanSTARRS/');
catfile=strcat(catdir,dt.name,'.txt');
M = csvread(catfile,1);
 
x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

sr = ((7./3600.0)*(pi/180.0)).^2;

    
if inst==1
    mps_arr = squeeze(M(:,16)');

else
    mps_arr = squeeze(M(:,17)');

end

sp=find(x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);
mps_arr = mps_arr(sp);
xps_arr = x_arr(sp);
yps_arr = y_arr(sp);

catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PSC/');
catfile=strcat(catdir,dt.name,'.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;
if inst==1
    mtm_arr = squeeze(M(:,8)');
else
    mtm_arr = squeeze(M(:,9)');
end
match_arr = squeeze(M(:,11)');
sp=find(match_arr==1 & ...
    x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);
mtm_arr = mtm_arr(sp);
xtm_arr = x_arr(sp);
ytm_arr = y_arr(sp);

%%% choose Ith s.t. masking fraction
Ith = 0.3;
mask = ones(1024);
mask(map > Ith) = 0;


%%% find radius
rtm_arr = zeros(size(mtm_arr));
for i = 1:numel(mtm_arr)
radmap=make_radius_map(mask,xtm_arr(i),ytm_arr(i));
r = 1;
while all(mask(radmap < r) == 0)
    r = r + 1;
end
rtm_arr(i) = r - 1;
end

rps_arr = zeros(size(mps_arr));
for i = 1:numel(mps_arr)
radmap=make_radius_map(mask,xps_arr(i),yps_arr(i));
r = 1;
while all(mask(radmap < r) == 0)
    r = r + 1;
end
rps_arr(i) = r - 1;

if rem(i, 1000)==0
    fprintf('find radius for PS source %d / %d\n',i,numel(mps_arr));
end
end

%%% 

dm = 1;
mbinsmin = round(min([mtm_arr,mps_arr]), 1) - dm;
mbinsmax = round(max([mtm_arr,mps_arr]), 1) + dm;
mbinedges = mbinsmin:dm:mbinsmax;
m_arr = (mbinedges(2:end) + mbinedges(1:end-1)) ./ 2;
r_arr = zeros(size(m_arr));
for i = 1:numel(m_arr)
    sp = find((mtm_arr) >= mbinedges(i) & (mtm_arr) < mbinedges(i+1));
    rtmi = rtm_arr(sp);
    sp = find((mps_arr) >= mbinedges(i) & (mps_arr) < mbinedges(i+1));
    rpsi = rps_arr(sp);
    if numel([rtmi, rpsi]) > 0
        r_arr(i) = min([rtmi, rpsi]);
    end

end

maskrad(ifield).maskfrac = numel(find(mask==1)) / 1024 / 1024;
maskrad(ifield).mtm_arr = mtm_arr;
maskrad(ifield).rtm_arr = rtm_arr;
maskrad(ifield).mps_arr = mps_arr;
maskrad(ifield).rps_arr = rps_arr;
maskrad(ifield).r_arr =r_arr;
maskrad(ifield).m_arr =m_arr;
end
%%
for ifield=4:8
    r_arr = maskrad(ifield).r_arr;
    m_arr = maskrad(ifield).m_arr;
    sp = (r_arr > 0) & (m_arr > 0) & (m_arr < 21);
    p = polyfit(m_arr(sp),r_arr(sp),4);
    fprintf('%.8e %.8e %.8e %.8e %.8e\n',p(1),p(2),p(3),p(4),p(5))
end 
%%% write the values to get_mask_radius.m
%%
for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    maskfrac = maskrad(ifield).maskfrac;
    mtm_arr = maskrad(ifield).mtm_arr;
    rtm_arr = maskrad(ifield).rtm_arr;
    mps_arr = maskrad(ifield).mps_arr;
    rps_arr = maskrad(ifield).rps_arr;
    r_arr = maskrad(ifield).r_arr;
    m_arr = maskrad(ifield).m_arr;
    
    sp = (r_arr > 0) & (m_arr > 0) & (m_arr < 21);

    rfit_arr = get_mask_radius(inst,ifield,m_arr);
    rfit_arr = rfit_arr./7;
    
    figure
    setwinsize(gcf,900,300)

    subplot(1,2,1)
    plot(mtm_arr, rtm_arr, '.', 'markersize',10);hold on
    plot(mps_arr, rps_arr, '.');
    plot(m_arr(sp), r_arr(sp),'ko', 'markersize',5);
    plot(m_arr, rfit_arr, 'k');
    xlabel('m_{AB}');
    ylabel('radius (pix)');
    legend({'2MASS', 'PanSTARRS', 'bin min', 'fit'})
    title(dt.name);
    subplot(1,2,2)
    semilogy(mtm_arr, rtm_arr, '.', 'markersize',10);hold on
    plot(mps_arr, rps_arr, '.');
    plot(m_arr(sp), r_arr(sp),'ko', 'markersize',5);
    plot(m_arr, rfit_arr, 'k');
    xlabel('m_{AB}');
    ylabel('radius (pix)');
end
%% making mask
flight = 40030;
inst = 2;
mypaths=get_paths(flight);

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

if inst == 1
    band = 'I';
else
    band = 'H';
end
[tmmask,tmnum] = make_mask_2m(flight,inst,band,ifield,1,13);
[psmask,psnum] = make_mask_ps(flight,inst,band,ifield,0,0,23);
strmask = psmask.*tmmask;
strnum = psnum + tmnum;

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');
maskdat.mask(ifield).strmask_stack = strmask;
maskdat.mask(ifield).strnum_stack = strnum;
save(strcat(loaddir,'maskdat'),'maskdat');
end
%%
for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
map = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps1_all.fits'));

figure
strmask1 = maskdat.mask(ifield).strmask;
imageclip(strmask1.*map);
title(sprintf('old mask %.1f %% unmasked', numel(find(strmask1))/1024/1024 * 100));
v = caxis;

figure
strmask = maskdat.mask(ifield).strmask_stack;
imageclip(strmask.*map);
title(sprintf('new mask %.1f %% unmasked', numel(find(strmask))/1024/1024 * 100));
caxis = v;

figure
strmask = maskdat.mask(ifield).strmask_stack;
histogram(strmask.*map);

end
%%
ifield = 5;
dt=get_dark_times(flight,inst,ifield);

srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
map = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps1_all.fits'));
strmask = maskdat.mask(ifield).strmask_stack;

figure
setwinsize(gcf,800,300)
subplot(1,2,2)
imageclip(map.*strmask);
title('elat30 sim map masked');
v = caxis;
c = colorbar;
c.Label.String = 'nW/m^2/sr';

subplot(1,2,1)
imageclip(map);
title('elat30 sim map');
caxis(v);
c = colorbar;
c.Label.String = 'nW/m^2/sr';
