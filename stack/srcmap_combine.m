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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_min_arr = [-5, 18, 19,   20, 20.5,   21, 21.5,   22, 22.5, 23, 24, 25];
m_max_arr = [18, 19, 20, 20.5,   21, 21.5,   22, 22.5,   23, 24, 25, 40];
group     = [ 1,  1,  1,    2,    2,    2,    2,    3,    3,  3,  3,  3];

map1 = zeros(720);
map2 = zeros(720);
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    if group(im) == 1
        map = fits_read(strcat(savedir,'unisphere_ihlmap',...
            num2str(m_min),'_',num2str(m_max),'sides1.fits'));
        map1 = map1 + map;
    elseif group(im)==2
        map = fits_read(strcat(savedir,'unisphere_ihlmap',...
            num2str(m_min),'_',num2str(m_max),'sides1.fits'));
        map2 = map2 + map;
    end
    
end
fits_write(strcat(savedir,'unisphere_ihlmap_sim1_all1'),map1);
fits_write(strcat(savedir,'unisphere_ihlmap_sim2_all1'),map2);

map1 = zeros(720);
map2 = zeros(720);
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    if group(im) == 1
        map = fits_read(strcat(savedir,'unisphere_ihlmap',...
            num2str(m_min),'_',num2str(m_max),'sides2.fits'));
        map1 = map1 + map;
    elseif group(im)==2
        map = fits_read(strcat(savedir,'unisphere_ihlmap',...
            num2str(m_min),'_',num2str(m_max),'sides2.fits'));
        map2 = map2 + map;
    end
    
end
fits_write(strcat(savedir,'unisphere_ihlmap_sim1_all2'),map1);
fits_write(strcat(savedir,'unisphere_ihlmap_sim2_all2'),map2);

map1 = zeros(720);
map2 = zeros(720);
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);

    if group(im) == 1
        map = fits_read(strcat(savedir,'unisphere_ihlmap',...
            num2str(m_min),'_',num2str(m_max),'sides3.fits'));
        map1 = map1 + map;
    elseif group(im)==2
        map = fits_read(strcat(savedir,'unisphere_ihlmap',...
            num2str(m_min),'_',num2str(m_max),'sides3.fits'));
        map2 = map2 + map;
    end
    
end
fits_write(strcat(savedir,'unisphere_ihlmap_sim1_all3'),map1);
fits_write(strcat(savedir,'unisphere_ihlmap_sim2_all3'),map2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_min_arr = [-5, 18, 19,   20, 20.5, 21, 22, 23, 24, 25];
m_max_arr = [18, 19, 20, 20.5,   21, 22, 23, 24, 25, 40];
group     = [ 1,  1,  1,    2,    2,  2,  2,  3,  3,  4];
for field_c={'UD_COSMOS','D_COSMOS','D_DEEP2-3','D_ELAISN1',...
        'XXM_00','XXM_01','XXM_10','XXM_11',...
        'XXM_20','XXM_21','XXM_30','XXM_31'}
    field=char(field_c);
    map1 = zeros(642);
    map2 = zeros(642);
    map3 = zeros(642);
    map4 = zeros(642);
    for im=1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
        if group(im) == 1
            map = fits_read(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'hsc_',field,'.fits'));
            map1 = map1 + map;
        elseif group(im)==2
            map = fits_read(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'hsc_',field,'.fits'));
            map2 = map2 + map;
        elseif group(im)==3
            map = fits_read(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'hsc_',field,'.fits'));
            map3 = map3 + map;
        elseif group(im)==4
            map = fits_read(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'hsc_',field,'.fits'));
            map4 = map4 + map;
        end

    end
    fits_write(strcat(savedir,'srcmap_sim1_all_hsc_',field),map1);
    fits_write(strcat(savedir,'srcmap_sim2_all_hsc_',field),map2);
    fits_write(strcat(savedir,'srcmap_sim3_all_hsc_',field),map3);
    fits_write(strcat(savedir,'srcmap_sim4_all_hsc_',field),map4);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_min_arr = [-5, 18, 19,   20, 20.5, 21, 22, 23, 24, 25];
m_max_arr = [18, 19, 20, 20.5,   21, 22, 23, 24, 25, 40];
group     = [ 1,  1,  1,    2,    2,  2,  2,  2,  2,  2];
for hsc_idx=0:11
    [field,~] = HSC_fields_info(hsc_idx);
    map1 = zeros(1024);
    map2 = zeros(1024);
    for im=1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
        if group(im) == 1
            map = fits_read(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'hsc_',field,'.fits'));
            map1 = map1 + map;
        elseif group(im)==2
            map = fits_read(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'hsc_',field,'.fits'));
            map2 = map2 + map;
        end

    end
    fits_write(strcat(savedir,'srcmap_sim1_all_hsc_',field),map1);
    fits_write(strcat(savedir,'srcmap_sim2_all_hsc_',field),map2);
end

return