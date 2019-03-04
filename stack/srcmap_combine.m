function srcmap_combine(flight, inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%run this on spire, and only download *all* in srcmap dir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mypaths=get_paths(flight);
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
savedir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_min_arr = [0,8:22];
m_max_arr = [8:23];
m_min_arr_2m = [0,7:21];
m_max_arr_2m = [7:22];
for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    psmap_raw = zeros(1024);
    for im=1:numel(m_min_arr)
        psmapsi = fits_read(strcat(srcmapdir,dt.name,'_srcmaps',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps.fits'));
        psmapgi = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps.fits'));
        psmapui = fits_read(strcat(srcmapdir,dt.name,'_srcmapu',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps.fits'));
        psmap_raw = psmap_raw + psmapsi + psmapgi + psmapui;
    end
    tmmap_raw = zeros(1024);
    for im=1:7
        m_min = m_min_arr_2m(im);
        m_max = m_max_arr_2m(im);
        tmmap = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'tm_PSsup.fits'));
        tmmap_raw = tmmap_raw + tmmap;
    end
    psmap_raw = psmap_raw + tmmap_raw;
    fits_write(strcat(savedir,dt.name,'_srcmap_ps_all'),psmap_raw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_min_arr = [0,8:22];
m_max_arr = [8:23];
m_min_arr_2m = [0,7:21];
m_max_arr_2m = [7:22];
for ifield = 4:8
    dt=get_dark_times(flight,inst,ifield);
    psmap_raw = zeros(1024);
    for im=1:numel(m_min_arr)
        psmapsi = fits_read(strcat(srcmapdir,dt.name,'_srcmaps',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps1.fits'));
        psmapgi = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps1.fits'));
        psmapui = fits_read(strcat(srcmapdir,dt.name,'_srcmapu',...
                num2str(m_min_arr(im)),'_',num2str(m_max_arr(im)),'ps1.fits'));
        psmap_raw = psmap_raw + psmapsi + psmapgi + psmapui;
    end
    tmmap_raw = zeros(1024);
    for im=1:7
        m_min = m_min_arr_2m(im);
        m_max = m_max_arr_2m(im);
        tmmap = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'tm_PSsup1.fits'));
        tmmap_raw = tmmap_raw + tmmap;
    end
    psmap_raw = psmap_raw + tmmap_raw;
    fits_write(strcat(savedir,dt.name,'_srcmap_ps1_all'),psmap_raw);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_min_arr = [-5, 18, 19,   20, 20.5, 21, 22, 23, 24, 25];
m_max_arr = [18, 19, 20, 20.5,   21, 22, 23, 24, 25, 40];
group     = [ 1,  1,  1,    2,    2,  2,  3,  3,  3,  3];

for ifield=4:8
    map1 = zeros(720);
    map2 = zeros(720);
    map3 = zeros(720);
    dt=get_dark_times(flight,inst,ifield);
    for im=1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
               
        maps = fits_read(strcat(savedir,dt.name,'_srcmap',...
            num2str(m_min),'_',num2str(m_max),'trilegal.fits'));
        mapg = fits_read(strcat(savedir,dt.name,'_srcmap',...
            num2str(m_min),'_',num2str(m_max),'sides.fits'));
        if group(im) == 1
            map1 = map1 + maps + mapg;
        elseif group(im)==2
            map2 = map2 + maps + mapg;
        else
            map3 = map3 + maps + mapg;
        end
    end
    fits_write(strcat(savedir,dt.name,'_srcmap_sim1_all'),map1);
    fits_write(strcat(savedir,dt.name,'_srcmap_sim2_all'),map2);
    fits_write(strcat(savedir,dt.name,'_srcmap_sim3_all'),map3);
end

return