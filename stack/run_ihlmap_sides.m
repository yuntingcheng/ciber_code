flight=40030;
inst=1;
interp = 1;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackmapdat',loaddir),'stackmapdat');
load(sprintf('%s/excessdat',loaddir),'excessdat');

savedir = strcat(mypaths.ciberdir, 'doc/20171018_stackihl/ihl_sides/TM',...
    num2str(inst),'/');

ifield=8;
dt=get_dark_times(flight,inst,ifield);
%% run all IHL
for m_min = 13:28
m_max = m_min + 1;
params = excessdat.all_avg.fit_params;
radius = excessdat.all_avg.fit_radius;

ihlmap = make_ihl_pl2_sides(flight,inst,ifield,interp,m_min,m_max,params,radius);
fits_write(strcat(savedir,dt.name,'_ihlmap_all',...
num2str(m_min),'_',num2str(m_max)),ihlmap);
end
%% run srcmap, ihlmap fraction of source 16<m<20
params = excessdat.all_avg.fit_params;
radius = excessdat.all_avg.fit_radius;

for im = 10:13
m_min = stackmapdat(ifield).m_min_arr(im);
m_max = stackmapdat(ifield).m_max_arr(im);
N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4; % CIBER field 4 deg^2
N_stack = stackmapdat(ifield).count_stackg_arr(im); % CIBER field 4 deg^2
frac = N_stack/N_helg;

[srcmap,ihlmap,mask] = make_partialmap_sides...
    (flight,inst,ifield,m_min,m_max,interp,params,radius,frac);
fits_write(strcat(savedir,dt.name,'_ihlmap_frac',...
num2str(m_min),'_',num2str(m_max)),ihlmap);
fits_write(strcat(savedir,dt.name,'_srcmap_frac',...
num2str(m_min),'_',num2str(m_max)),srcmap);
fits_write(strcat(savedir,dt.name,'_mask_frac',...
num2str(m_min),'_',num2str(m_max)),mask);
end
