function run_ihlmap_sim(flight, inst, rvir, im)
mypaths=get_paths(flight);

m_min_arr = [-5, 18, 19,   20, 20.5, 21, 21.5, 22, 22.5, 23, 24, 25];
m_max_arr = [18, 19, 20, 20.5,   21, 21.5, 22, 22.5, 23, 24, 25, 40];

savedir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

m_min = m_min_arr(im);
m_max = m_max_arr(im);

map = make_ihlmap_unisphere(flight,inst,1,m_min,m_max,rvir);
fits_write(strcat(savedir,'unisphere_ihlmap',...
    num2str(m_min),'_',num2str(m_max),'sides',num2str(rvir)),map);

return