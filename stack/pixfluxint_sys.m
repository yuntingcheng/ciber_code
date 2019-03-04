%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%see how the systematics in PSF bias the ratio
%i.e. use wrong PSF or wrong ra/dec for optimal photometry
% systematics: 
%1. 1.2x wider PSF
%2. 0.5 pix offset of ra dec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
field=4;
quad_arr=['A'];
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/simflux/TM',...
    num2str(inst),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/simflux/plots/TM',...
    num2str(inst),'/');

dt=get_dark_times(flight,inst,field);
ukcatlab='UKIDSS input [Jy]';
uk0simlab='UKIDSS sim true [Jy]';
uk1simlab='UKIDSS sim 1.2$\sigma$ [Jy]';
uk2simlab='UKIDSS sim 0.5 pix offset [Jy]';
cbdatlab='CIBER data [Jy]';
%% convert I to flux
sr=((7./3600.0)*(pi/180.0))^2;
if inst==1
    wleff=1.05e-6;% [m]
else
    wleff=1.79e-6;% [m]
end
nueff=3e8/wleff;% [Hz]
I2F=sr/nueff*1e-9*1e26; %[Jy]
%%
npix=3;
for quad=quad_arr
figure
setwinsize(gcf,1000,1000)
for i=1:4
subplot(2,2,i)
fname=strcat(loaddir,dt.name,'_',quad,'_pix',num2str(npix),'_flux_sys.txt');
M = csvread(fname);
Mcat_arr=M(:,6)';
Icat_arr=M(:,7)';
Iuk0_arr=M(:,8)';
Iuk1_arr=M(:,9)';
Iuk2_arr=M(:,10)';

fname=strcat(loaddir,dt.name,'_',quad,'_pix',num2str(npix),'_flux_opt.txt');
M = csvread(fname);
Icb_arr=M(:,8)';

if i==1
    Ix_arr=Icat_arr;xlab=ukcatlab;
    Iy_arr=Iuk0_arr;ylab=uk0simlab;
    sname='uk0_in';
elseif i==2
    Ix_arr=Icat_arr;xlab=ukcatlab;
    Iy_arr=Iuk1_arr;ylab=uk1simlab;
    sname='uk1_in';
elseif i==3
    Ix_arr=Icat_arr;xlab=ukcatlab;
    Iy_arr=Iuk2_arr;ylab=uk2simlab;
    sname='uk2_in';
else
    Ix_arr=Icat_arr;xlab=ukcatlab;
    Iy_arr=Icb_arr;ylab=cbdatlab;
    sname='uk2_in';
    
end

Iedge_arr=logspace(log10(30),log10(1e6),21);
Iedge_arr=[Iedge_arr,1e20];

x_arr=zeros(1,numel(Iedge_arr)-1);
ex_arr=zeros(1,numel(Iedge_arr)-1);
y_arr=zeros(1,numel(Iedge_arr)-1);
ey_arr=zeros(1,numel(Iedge_arr)-1);
for iI=1:numel(Iedge_arr)-1
    sp=find(Ix_arr>Iedge_arr(iI) & Ix_arr<=Iedge_arr(iI+1));
    yuse=Iy_arr(sp);
    x_arr(iI)=nanmedian(Ix_arr(sp));
    ex_arr(iI)=std(Ix_arr(sp))./sqrt(numel(sp));
    y_arr(iI)=nanmedian(yuse);
    ey_arr(iI)=std(yuse)./sqrt(numel(sp));
    
end
x_arr=x_arr.*I2F;
y_arr=y_arr.*I2F;
ex_arr=ex_arr.*I2F;
ey_arr=ey_arr.*I2F;

errorbar(x_arr,y_arr,ey_arr,ey_arr,ex_arr,ex_arr,'.');hold on

set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([1-6,1e0])
ylim([1e-6,1e0])
plot([1e-6,1e0],[1e-6,1e0],'k');
title(strcat(dt.name,'\_',quad,'\_',num2str(npix),'x',num2str(npix)),'fontsize',15)
xlabel(xlab,'interpreter','latex','fontsize',18)
ylabel(ylab,'interpreter','latex','fontsize',18)

end
imname=strcat(savedir,dt.name,'_sys_',quad);
print(imname,'-dpng');%close
end

