flight=40030;
inst=2;
ifield=8;
%m_min = 13;
m_max = m_min + 1;

mypaths=get_paths(flight);

savedir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
interp = 1;
dt=get_dark_times(flight,inst,ifield);
maps = make_srcmap_sides(flight,inst,ifield,m_min,m_max,interp);
fits_write(strcat(savedir,dt.name,'_srcmap',...
    num2str(m_min),'_',num2str(m_max),'sides'),maps);
%% plot the source map GIF
flight=40030;
inst=2;
ifield=8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

pltsavedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));
gifname=strcat(pltsavedir,dt.name,'_srcmap.gif');

srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

dt=get_dark_times(flight,inst,ifield);
maptot = zeros(720);
m_arr = [];
mean_arr = [];

fig = figure;
m_min_arr = 13:28;
for im = 1:numel(m_min_arr)
    m_min = m_min_arr(im);
    m_max = m_min + 1;
    m_arr  = [m_arr mean([m_min,m_max])];
    map = fits_read(strcat(srcmapdir,dt.name,'_srcmap',...
        num2str(m_min),'_',num2str(m_max),'sides.fits'));
    maptot = map + maptot;

    mean_arr = [mean_arr, mean(map(:))];
    imageclip(log10(maptot));
    h = colorbar;
    caxis([-0.25,1.7]);
    title(sprintf('m_{AB}^y < %d',m_max),'fontsize',20);
    ylabel(h,'log10(nW/m^2/sr)','fontsize',20);
    
    % make gif
    frame = getframe(fig);
    [A,map] = rgb2ind(frame2im(frame),256);
    if im == 1
        imwrite(A,map,gifname,'gif','LoopCount',Inf,'DelayTime',0.5);
    else
        imwrite(A,map,gifname,'gif','WriteMode','append','DelayTime',0.5);
    end

end

semilogy(m_arr+0.5,mean_arr,'o-');hold on
semilogy(m_arr,cumsum(mean_arr),'o-');
ylim([0.08,20])
xlabel('m_{AB}^y','fontsize',20);
ylabel('<I>(nW/m^2/sr)','fontsize',20);
legend({'d<I>/dm','<I>(<m)'},'fontsize',20,'location','northwest')
savename=strcat(pltsavedir,dt.name,'_meanI');
print(savename,'-dpng');%close