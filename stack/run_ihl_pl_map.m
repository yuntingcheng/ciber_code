flight=40030;
inst=1;
ifield=8;
dt=get_dark_times(flight,inst,ifield);
srcmapdir = strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
loaddir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/');
load(sprintf('%s/ciber_ps/%s_profstack',loaddir,dt.name));
savedir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_pl/TM',...
    num2str(inst),'/');

m_min_arr = [0,8:22];
m_max_arr = [8:23];
%% run the IHL maps 
%{
for im=10:12
    m_min = profstack(im).m_min;
    m_max = profstack(im).m_max;
    plparams = profstack(im).pl_fit_params;
    ihlmap = make_ihl_pl(flight,inst,ifield,m_min,m_max,plparams);
    
    fits_write(strcat(savedir,dt.name,'_ihlmap_pl',...
    num2str(m_min),'_',num2str(m_max)),ihlmap);

end
%}
%% 
pixscale=7;
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

for im=10:12
    
    m_min = profstack(im).m_min;
    m_max = profstack(im).m_max;
    
    
    psmapg = fits_read(strcat(srcmapdir,dt.name,'_srcmapg',...
         num2str(m_min),'_',num2str(m_max),'ps.fits'));
    ihlmap = fits_read(strcat(savedir,dt.name,'_ihlmap_pl',...
    num2str(m_min),'_',num2str(m_max),'.fits'));
    ihlmap(find(ihlmap==Inf)) = 0;
    
    mask = make_mask_ps(flight,inst,ifield,1,m_min,m_max);
    mkk =get_mkk_sim(mask,pixscale,binl,10,numel(binl),1,ones(1024),0,NaN);

    mmap = (psmapg+ihlmap).*mask;
    [Clihl,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));
    
    mmap = (psmapg).*mask;
    [Clpsf,~,~,l]=get_Cl(mmap,mask,mkk,pixscale,ones(1024));

    loglog(l,l.*(l+1).*Clpsf./2./pi,'--','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF'));hold on
    loglog(l,l.*(l+1).*Clihl./2./pi,'-','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m_min),'<m<',num2str(m_max),' PSF+IHL'));
    
end
xlim([1e2,2e5]);
ylim([1e-6,1e-2]);
title(sprintf('power-law IHL'),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_pl/TM',...
    num2str(inst),'/'));
savename=strcat(savedir,dt.name,'_ps');
print(savename,'-dpng');%close
