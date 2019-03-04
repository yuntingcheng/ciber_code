flight = 40030;
inst = 1;
ifield = 5;
%%
mypaths=get_paths(flight);

m_min_arr = [0,8:22];
m_max_arr = [8:23];

m_min_arr_2m = [0,7:21];
m_max_arr_2m = [7:22];

dt=get_dark_times(flight,inst,ifield);

%%% get srcmap %%%
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
psmap = zeros(1024);
for im=1:numel(m_min_arr)
    psmapsi = fits_read(strcat(srcmapdir,dt.name,'_srcmaps',...
            num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps1.fits'));
    psmapgi = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
            num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps1.fits'));
    psmapui = fits_read(strcat(srcmapdir,dt.name,'_srcmapu',...
            num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps1.fits'));
    psmap = psmap + psmapsi + psmapgi + psmapui;
end
tmmap = zeros(1024);
for im=1:7
    m_min = m_min_arr_2m(im);
    m_max = m_max_arr_2m(im);
    tmmapi = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
            num2str(m_min),'_',num2str(m_max),'tm_PSsup1.fits'));
    tmmap = tmmap + tmmapi;
end
map = psmap + tmmap;

%%% get cats %%%
loaddir = '/Users/ytcheng/Downloads/';
catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PanSTARRS/');
catfile=strcat(catdir,dt.name,'.txt');
M = csvread(catfile,1);
 
x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;
 
my_arr=squeeze(M(:,9)');
cls_arr=squeeze(M(:,10)');
 
sr = ((7./3600.0)*(pi/180.0)).^2;

    
if inst==1
    mlin_arr = squeeze(M(:,11)');
else
    mlin_arr = squeeze(M(:,12)');
end
[mps_arr, ~] = get_corrected_mag(inst, mlin_arr, my_arr, cls_arr);

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
%% choose Ith s.t. masking fraction ~ 0.35
Ith = 3;
mask = ones(1024);
mask(map > Ith) = 0;
numel(find(mask==1)) / 1024 / 1024
%% find radius
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
%% 
dm = 0.2;
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

sp = (r_arr > 0) | (m_arr > 20);
p = polyfit(m_arr(sp),r_arr(sp),4);

% copy this poly params to get_mask_radius.m
fprintf('%.8e %.8e %.8e %.8e %.8e\n',p(1),p(2),p(3),p(4),p(5))
%%
rfit_arr = get_mask_radius(inst,ifield,m_arr);
rfit_arr = rfit_arr./7;
%%
figure
setwinsize(gcf,900,300)

subplot(1,2,1)
plot(mtm_arr, rtm_arr, '.', 'markersize',10);hold on
plot(mps_arr, rps_arr, '.');
plot(m_arr(sp), r_arr(sp),'ko', 'markersize',10);
plot(m_arr, rfit_arr, 'k');
legend({'2MASS', 'PanSTARRS', 'bin min', 'fit'})

subplot(1,2,2)
semilogy(mtm_arr, rtm_arr, '.', 'markersize',10);hold on
plot(mps_arr, rps_arr, '.');
plot(m_arr(sp), r_arr(sp),'ko', 'markersize',10);
plot(m_arr, rfit_arr, 'k');
legend({'2MASS', 'PanSTARRS', 'bin min', 'fit'})
%%
if inst == 1
    band = 'I';
else
    band = 'H';
end
[tmmask,tmnum] = make_mask_2m(flight,inst,band,ifield,1,13);
[psmask,psnum] = make_mask_ps(flight,inst,band,ifield,0,0,23);
strmask = psmask.*tmmask;
strnum = psnum + tmnum;
%%
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');

figure
strmask1 = maskdat.mask(ifield).strmask;
imageclip(strmask1.*map);
title(sprintf('old mask %.1f %% unmasked', numel(find(strmask1))/1024/1024 * 100));
v = caxis;
figure
imageclip(strmask.*map);
caxis = v;
title(sprintf('new mask %.1f %% unmasked', numel(find(strmask))/1024/1024 * 100));
%%
maskdat.mask(ifield).strmask_stack = strmask;
maskdat.mask(ifield).strnum_stack = strnum;

% try the aggressive strmask (~ 20% unmasked)
% maskdat.mask(ifield).strmask_stack_aggressive = strmask;
% maskdat.mask(ifield).strnum_stack_aggressive = strnum;

save(strcat(loaddir,'maskdat'),'maskdat');