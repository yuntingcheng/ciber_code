function write_ihl_sim_coord(flight,inst,ifield)
mypaths=get_paths(flight);

savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');

dt=get_dark_times(flight,inst,ifield);

countg_arr = stackmapdat(ifield).count_stackg_arr;
x_arr = [];
y_arr = [];
xr_arr = [];
yr_arr = [];
m_arr = [];
use_arr = [];

[~,l,lbin,Clin,~]=sigCl_extrap_mz14;
sigmap = map_from_power_spec(lbin,Clin,1024,1024,7,1);
sigmap=sigmap./std(sigmap(:));
sigmap(find(sigmap<-1))=-1;
sigmap=sigmap+1;

for im=10:13
    m_min = stackmapdat(ifield).m_min_arr(im);
    m_max = stackmapdat(ifield).m_max_arr(im);

    N_helg = IGLcounts_helgason(inst,(m_min+m_max)/2) * 4;
    N_helg = round(N_helg);
    N_stack = countg_arr(im);


    sp = find(sigmap);
    ind=randsample(sp,N_helg,true,sigmap(sp));
    [x,y] = ind2sub(size(sigmap),ind);
    x = x + rand(size(x)) - 0.5;
    y = y + rand(size(y)) - 0.5;
    
    xr = rand(size(x))*1024 + 0.5;
    yr = rand(size(y))*1024 + 0.5;
    
    m = rand(size(y)) + m_min;
    
    ind=randsample(1:N_helg,N_stack);
    use = zeros(size(x));
    use(ind) = 1;
    
    x_arr = [x_arr x'];
    y_arr = [y_arr y'];
    xr_arr = [xr_arr xr'];
    yr_arr = [yr_arr yr'];
    m_arr = [m_arr m'];
    use_arr = [use_arr use'];    
end

M = zeros(numel(x_arr),6);
M(:,1) = x_arr;
M(:,2) = y_arr;
M(:,3) = xr_arr;
M(:,4) = yr_arr;
M(:,5) = m_arr;
M(:,6) = use_arr;

csvwrite(strcat(savedir,dt.name,'_simcoord.csv'),M);
return