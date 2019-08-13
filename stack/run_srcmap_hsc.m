function run_srcmap_hsc(flight, inst, hsc_idx)
mypaths=get_paths(flight);

m_min_arr = [-5, 18, 19,   20, 20.5, 21, 22, 23, 24, 25];
m_max_arr = [18, 19, 20, 20.5,   21, 22, 23, 24, 25, 40];

usePSF_arr= (m_min_arr < 20);

savedir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

ifield = 8;
dt=get_dark_times(flight,inst,ifield);
[field,~] = HSC_fields_info(hsc_idx);
for im=1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    
    map = make_srcmap_hsc(flight,inst,ifield,field,m_min,m_max,usePSF_arr(im));
    
    fits_write(strcat(savedir,dt.name,'_srcmap',...
        num2str(m_min),'_',num2str(m_max),'hsc_',field),map);
end

return