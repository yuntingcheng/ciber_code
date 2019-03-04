%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the total source flux of CIBER, 2MASS sim map, 2M catalog flux. 
% Using both optimal and aperture photometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
quad_arr=['A','B','C','D'];
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/simflux/TM',...
    num2str(inst),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/simflux/plots/TM',...
    num2str(inst),'/');

ukcatlab='2MASS input [nW/m2/sr]';
cbdatlab='CIBER data [nW/m2/sr]';

%%% convert I to flux
sr=((7./3600.0)*(pi/180.0))^2;
if inst==1
    wleff=1.05e-6;% [m]
    rcal=0.58;
else
    wleff=1.79e-6;% [m]
    rcal=0.38;
end
nueff=3e8/wleff;% [Hz]
I2F=sr/nueff*1e-9*1e26; %[Jy]
%% plot 3x3 different data combine quad
for ifield=4%4:8
dt=get_dark_times(flight,inst,ifield);

npix=3;

Mcat_arr=[];
Icat_arr=[];
Icb_arr=[];

for quad=quad_arr
fname=strcat(loaddir,dt.name,'_',quad,'_pix',num2str(npix),'_flux_opt_2m.txt');
M = csvread(fname);

Mcat_arr=[Mcat_arr M(:,6)'];
Icat_arr=[Icat_arr M(:,7)'];
Icb_arr=[Icb_arr M(:,10)'];
end


Ix_arr=Icat_arr;xlab=ukcatlab;
Iy_arr=Icb_arr;ylab=cbdatlab;
sname='cb_in';


Iedge_arr=linspace(0,2000,40);

x_arr=zeros(1,numel(Iedge_arr)-1);
ex_arr=zeros(1,numel(Iedge_arr)-1);
y_arr=zeros(1,numel(Iedge_arr)-1);
ey_arr=zeros(1,numel(Iedge_arr)-1);
count_arr=zeros(1,numel(Iedge_arr)-1);
for iI=1:numel(Iedge_arr)-1
    sp=find(Ix_arr>Iedge_arr(iI) & Ix_arr<=Iedge_arr(iI+1));
    yuse=Iy_arr(sp);
    x_arr(iI)=nanmedian(Ix_arr(sp));
    ex_arr(iI)=std(Ix_arr(sp))./sqrt(numel(sp));
    y_arr(iI)=nanmedian(yuse);
    ey_arr(iI)=std(yuse)./sqrt(numel(sp));
    count_arr(iI)=numel(sp);
    
end

errorbar(x_arr,y_arr./rcal,ey_arr./rcal,ey_arr./rcal,ex_arr,ex_arr,'r.');hold on
errorbar(x_arr,y_arr,ey_arr,ey_arr,ex_arr,ex_arr,'b.');
%set(gca,'xscale','log')
%set(gca,'yscale','log')
xlim([Iedge_arr(1),Iedge_arr(end)])
ylim([Iedge_arr(1),Iedge_arr(end)])
plot([Iedge_arr(1),Iedge_arr(end)],[Iedge_arr(1),Iedge_arr(end)],'k');
title(strcat(dt.name,'\_',num2str(npix),'x',num2str(npix)),'fontsize',15)
xlabel(xlab,'interpreter','latex','fontsize',18)
ylabel(ylab,'interpreter','latex','fontsize',18)

end
%imname=strcat(savedir,dt.name,'_data_',sname);
%print(imname,'-dpng');%close
