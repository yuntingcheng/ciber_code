%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparing the total source flux of CIBER, 2MASS sim map, 2M catalog flux. 
% Using both optimal and aperture photometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sum=1;% do summation not optimal photometry
flight=40030;
inst=1;
quad_arr=['A','B','C','D'];
npix_arr=[1,3,5,7];
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/simflux/TM',...
    num2str(inst),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/simflux/plots/TM',...
    num2str(inst),'/');

ukcatlab='2MASS input [Jy]';
uksimlab='2MASS sim [Jy]';
cbdatlab='CIBER data [Jy]';

%%% convert I to flux
sr=((7./3600.0)*(pi/180.0))^2;
if inst==1
    wleff=1.05e-6;% [m]
else
    wleff=1.79e-6;% [m]
end
nueff=3e8/wleff;% [Hz]
I2F=sr/nueff*1e-9*1e26; %[Jy]
%% plot 3x3 different data combine quad
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);

npix=3;

Mcat_arr=[];
Icat_arr=[];
Icb_arr=[];
Iuk_arr=[];
Icb1_arr=[];
Iuk1_arr=[];

for quad=quad_arr
fname=strcat(loaddir,dt.name,'_',quad,'_pix',num2str(npix),'_flux_opt_2m.txt');
M = csvread(fname);

Mcat_arr=[Mcat_arr M(:,6)'];
Icat_arr=[Icat_arr M(:,7)'];
Icb_arr=[Icb_arr M(:,8)'];
Iuk_arr=[Iuk_arr M(:,9)'];
Icb1_arr=[Icb1_arr M(:,10)'];
Iuk1_arr=[Iuk1_arr M(:,11)'];
end

if sum==1
    Icb_arr=Icb1_arr;Iuk_arr=Iuk1_arr;
end

figure
setwinsize(gcf,1000,1000)
for i=1:4
subplot(2,2,i)

if i==1
    Ix_arr=Icat_arr;Ix1_arr=Icat_arr;xlab=ukcatlab;
    Iy_arr=Iuk_arr;Iy1_arr=Iuk1_arr;ylab=uksimlab;
    sname='2m_in';

elseif i==2
    Ix_arr=Icat_arr;Ix1_arr=Icat_arr;xlab=ukcatlab;
    Iy_arr=Icb_arr;Iy1_arr=Icb1_arr;ylab=cbdatlab;
    sname='cb_in';

elseif i==3
    Ix_arr=Iuk_arr;Ix1_arr=Iuk1_arr;xlab=uksimlab;
    y_arr=Icb_arr;Iy1_arr=Icb1_arr;ylab=cbdatlab;
    name='cb_2m';
else
    Ix_arr=Icb_arr;xlab=cbdatlab;
    Iy_arr=Iuk_arr;ylab=uksimlab;
    sname='2m_cb';
end

Iedge_arr=logspace(log10(300),log10(1e6),21);
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
xlim([1e-4,1e1])
ylim([1e-4,1e1])
plot([1e-4,1e1],[1e-4,1e1],'k');
title(strcat(dt.name,'\_',num2str(npix),'x',num2str(npix)),'fontsize',15)
xlabel(xlab,'interpreter','latex','fontsize',18)
ylabel(ylab,'interpreter','latex','fontsize',18)

end
imname=strcat(savedir,dt.name,'_data_',sname);
if sum==1
    imname=strcat(savedir,dt.name,'_data_',sname,'_sum');
end
print(imname,'-dpng');%close
end
%% fit the scaling factor (do two fields in the same plots)
npix=3;
if inst==1
    r_arr=0.5:0.001:0.8;
else
    r_arr=0.3:0.001:0.6;
end
figure
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);

Mcat_arr=[];
Icat_arr=[];
Icb_arr=[];
Iuk_arr=[];
Icb1_arr=[];
Iuk1_arr=[];
for quad=quad_arr
fname=strcat(loaddir,dt.name,'_',quad,'_pix',num2str(npix),'_flux_opt_2m.txt');
M = csvread(fname);

Mcat_arr=[Mcat_arr M(:,6)'];
Icat_arr=[Icat_arr M(:,7)'];
Icb_arr=[Icb_arr M(:,8)'];
Iuk_arr=[Iuk_arr M(:,9)'];
Icb1_arr=[Icb1_arr M(:,10)'];
Iuk1_arr=[Iuk1_arr M(:,11)'];
end
if sum==1
    Icb_arr=Icb1_arr;Iuk_arr=Iuk1_arr;
end

Ix_arr=Icat_arr;Ix1_arr=Icat_arr;xlab=ukcatlab;
Iy_arr=Icb_arr;Iy1_arr=Icb1_arr;ylab=cbdatlab;
sname=strcat(dt.name);


Iedge_arr=logspace(log10(1e3),log10(1e5),21);
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

chi2_arr=zeros(size(r_arr));
for ir=1:numel(r_arr)
    r=r_arr(ir);
    chi2in=(y_arr-r.*x_arr).^2./(ey_arr.^2+r^2.*ex_arr.^2);
    chi2_arr(ir)=nansum(chi2in);
end
dof=numel(find(chi2in==chi2in))-1;
bestr=r_arr(find(chi2_arr==min(chi2_arr)));

plot(r_arr,chi2_arr./dof,'color',get_color(ifield-3),'linestyle','-',...
    'DisplayName',sprintf('%s, F_{CIBER}/F_{2M}=%.3f',sname,bestr));hold on


end
ylim([0,5])
xlim([r_arr(1), r_arr(end)])
xlabel('$F_{CIBER}/F_{2M}$','interpreter','latex','fontsize',18)
ylabel('$\chi^2/dof$','interpreter','latex','fontsize',18)
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
imname=strcat(savedir,'chi2_2m');
if sum==1
    imname=strcat(savedir,'chi2_2m_sum');
end

print(imname,'-dpng');%close

