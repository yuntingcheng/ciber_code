function run_srcmap(flight, inst, ifield_arr)
mypaths=get_paths(flight);

m_min_arr = [0,8:22];
m_max_arr = [8:23];

savedir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
% linear interp to CIBER band
% interp = 1;
% interp1 after correction & select with I/H band
interp = 2;
%%

for ifield=ifield_arr
    for im=1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
    
        dt=get_dark_times(flight,inst,ifield);
        maps = make_srcmap_ps(flight,inst,ifield,-1,m_min,m_max,interp);
        mapg = make_srcmap_ps(flight,inst,ifield,1,m_min,m_max,interp);
        mapu = make_srcmap_ps(flight,inst,ifield,2,m_min,m_max,interp);

        if interp == 1
            fits_write(strcat(savedir,dt.name,'_srcmaps',...
                num2str(m_min),'_',num2str(m_max),'ps'),maps);     
            fits_write(strcat(savedir,dt.name,'_srcmapg',...
                num2str(m_min),'_',num2str(m_max),'ps'),mapg);
            fits_write(strcat(savedir,dt.name,'_srcmapu',...
                num2str(m_min),'_',num2str(m_max),'ps'),mapu);
        elseif interp == 2
            fits_write(strcat(savedir,dt.name,'_srcmaps',...
                num2str(m_min),'_',num2str(m_max),'ps1'),maps);     
            fits_write(strcat(savedir,dt.name,'_srcmapg',...
                num2str(m_min),'_',num2str(m_max),'ps1'),mapg);
            fits_write(strcat(savedir,dt.name,'_srcmapu',...
                num2str(m_min),'_',num2str(m_max),'ps1'),mapu);
        end
    end
end
%%

%%%%%%%%%% 2MASS %%%%%%%%%%%%%%%%
m_min_arr = [0,7:21];
m_max_arr = [7:22];
PSmatch = 0;
for ifield=ifield_arr
    for im=1:numel(m_min_arr)
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
    
        dt=get_dark_times(flight,inst,ifield);
        map = make_srcmap_2m(flight,inst,ifield,m_min,m_max,interp,PSmatch);
        
        if interp == 1
            fits_write(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'tm'),map);
        elseif interp == 2
            fits_write(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'tm1'),map);
        end
    end
end

%%

m_min_arr = [0,7:21];
m_max_arr = [7:22];
PSmatch = 1;
for ifield=ifield_arr
    for im=1:7
        m_min = m_min_arr(im);
        m_max = m_max_arr(im);
    
        dt=get_dark_times(flight,inst,ifield);
        map = make_srcmap_2m(flight,inst,ifield,m_min,m_max,interp,PSmatch);
        if interp == 1
            fits_write(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'tm_PSsup'),map);
        elseif interp == 2
            fits_write(strcat(savedir,dt.name,'_srcmap',...
                num2str(m_min),'_',num2str(m_max),'tm_PSsup1'),map);
        end

    end
end

return