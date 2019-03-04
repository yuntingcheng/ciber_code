flight=40030;
inst=1;
ifield=8;
npix=800;
pixsize=7;

dt=get_dark_times(flight,inst,ifield);
srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');

savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_const/TM',...
    num2str(inst),'/'));
quad_arr=['A','B','C','D'];
m_arr=13:17;
load(strcat(savedir,'offsetdat'),'offsetdat');
rad=500;
%%

totsrcmap=zeros(1024);
totihlmap=zeros(1024);
totmask=ones(1024);
figure
for im=1:numel(m_arr)
    m=m_arr(im);
    amp=offsetdat(im).off;
    if im==4
        amp=(offsetdat(3).off+offsetdat(5).off)./2;
    end
    
    ihlmap=make_ihl_const(flight,inst,ifield,m,m+1,amp,rad);
    mask=make_galmask_2m(flight,inst,ifield,m,m+1);
    
    srcmapg=zeros(1024);
    for iquad=1:4
        quad=quad_arr(iquad);
        stmmapg=fits_read(strcat(srcmapdir,...
            dt.name,'_',quad,'_srcmapg',num2str(m),'_2m.fits'));

        if iquad==1
            srcmapg(1:512,1:512)=stmmapg;
        elseif iquad==2
            srcmapg(513:1024,1:512)=stmmapg;
        elseif iquad==3
            srcmapg(1:512,513:1024)=stmmapg;
        else
            srcmapg(513:1024,513:1024)=stmmapg;
        end
    end
    totsrcmap=totsrcmap+srcmapg;
    totihlmap=totihlmap+ihlmap;
    totmask=totmask.*mask;
    
    mkk=numel(find(mask))/1024^2;
    mmaps=(srcmapg-mean(srcmapg(find(mask)))).*mask;
    [Cls,l]=get_angular_spec(mmaps,mmaps,pixsize);
    Cls=Cls./mkk;
    srcmapgh=srcmapg+ihlmap;
    mmaph=(srcmapgh-mean(srcmapgh(find(mask)))).*mask;
    [Clh,l]=get_angular_spec(mmaph,mmaph,pixsize);
    Clh=Clh./mkk;


    loglog(l,l.*(l+1).*Cls./2./pi,'--','color',get_color(im),'linewidth',2,...
        'DisplayName',strcat(num2str(m),'<m<',num2str(m+1),' PSF'));hold on
    loglog(l,l.*(l+1).*Clh./2./pi,'-','color',get_color(im),'linewidth',2,...
                'DisplayName',strcat(num2str(m),'<m<',num2str(m+1),' PSF+IHL'));
    drawnow
    
end

mkk=numel(find(totmask))/1024^2;
mmaps=(totsrcmap-mean(totsrcmap(find(totmask)))).*totmask;
[Cls,l]=get_angular_spec(mmaps,mmaps,pixsize);
Cls=Cls./mkk;
totsrcmaph=totsrcmap+totihlmap;
mmaph=(totsrcmaph-mean(totsrcmaph(find(totmask)))).*totmask;
[Clh,l]=get_angular_spec(mmaph,mmaph,pixsize);
Clh=Clh./mkk;

loglog(l,l.*(l+1).*Cls./2./pi,'k--','linewidth',3,...
    'DisplayName','total PSF');hold on
loglog(l,l.*(l+1).*Clh./2./pi,'k-','linewidth',3,...
            'DisplayName','total PSF+IHL');

xlim([1e2,2e5]);
ylim([1e-5,1e-1]);
title(sprintf('radius=%d arcsec',rad),'fontsize',20);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northeastoutside');
set(h,'fontsize',12)

savename=strcat(savedir,dt.name,'_r',num2str(rad),'_PS');
print(savename,'-dpng');%close

%%%%%%%%%%%%%%%%%%%
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(totsrcmap);
title('gals*PSF');
v=caxis;
subplot(2,2,2)
imageclip(totihlmap+totsrcmap);
title('gals*PSF+IHL');
caxis(v);
subplot(2,2,3)
imageclip(totsrcmap.*totmask);
title('gals*PSF*mask');
caxis(v);
subplot(2,2,4)
imageclip((totihlmap+totsrcmap).*totmask);
title('(gals*PSF+IHL)*mask');
caxis(v);

set(gcf,'NextPlot','add');
axes; 
set(gca,'Visible','off'); 
h = title(sprintf('radius=%d arcsec',rad),'fontsize',20);
set(h,'Visible','on');

savename=strcat(savedir,dt.name,'_r',num2str(rad),'_map');
print(savename,'-dpng');%close
