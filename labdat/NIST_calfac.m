band=2;
nistdir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/NIST/TM'...
                    ,num2str(band),'/');                
maskdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskdir,'band',num2str(band),'_mask_inst'),'mask_inst');


ffdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(ffdir,'band',num2str(band),'_FF5'),'FF5');
framedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/NIST/TM'...
                    ,num2str(band),'/framedat/'); 
savedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/NIST/TM'...
                    ,num2str(band),'/plots/');                

switch band
case 1
nistlabel=[1,1,1,1,2,2,2,2,2,2];%different light level
time1={'10-40-20';'10-41-15'};
time2={'10-42-59';'10-43-23';'10-44-52';'10-45-13';'10-45-40'};
time3={'10-50-45';'10-51-10';'10-51-38';'10-52-12';'10-52-37'};
time4={'10-54-11';'10-55-10'};
time_arr={time1,time2,time3,time4};
% central freq from Zemcov+14 table S1
% range from Bock+13 table 2
wlmin=(0.9+(1.05-1.12))*1e3;
wlmax=(1.32+(1.05-1.12))*1e3;
case 2
nistlabel=[1,1,1,1,1,1,2,2,2,2,2,2];
time1={'11-03-27';'11-04-44';'11-05-04';'11-05-27';'11-05-54'};
time2={'11-08-36';'11-08-59';'11-09-18';'11-09-35';'11-09-56'};
time3={'11-13-30';'11-13-51';'11-14-10';'11-14-44';'11-15-06'};
time4={'11-16-35';'11-16-59';'11-17-19';'11-17-40';'11-18-03'};
time_arr={time1,time2,time3,time4};
% central freq from Zemcov+14 table S1
% range from Bock+13 table 2
wlmin=(1.15+(1.79-1.56))*1e3;
wlmax=(2.04+(1.79-1.56))*1e3;
wlmax=2000;%above 2000 are garbage
end
cp=get_cal_params('flight',40030);
cal=cp(band).apf2eps*cp(band).eps2nWpm2ps;
%% get NIST data
load(strcat(nistdir,'NISTdat'),'NISTdat');
nist_arr=zeros(1,numel(NISTdat));
figure
for i=1:numel(NISTdat)
Wavelength_arr=NISTdat(i).Wavelength_arr;%nm
L_arr=NISTdat(i).L_arr;%W/cm^3/sr
lIl_arr=(Wavelength_arr.*1e-7).*(L_arr.*1e4).*1e9;%nW/m^2/sr
nist_arr(i)=mean(lIl_arr(find(Wavelength_arr>wlmin & Wavelength_arr<wlmax)));

if nistlabel(i)==1;col='b';else col='r';end
plot(Wavelength_arr,lIl_arr,'color',col);hold on
xlim([600 2100])
end
xlabel('$\lambda(nm)$','interpreter','latex','fontsize',18)
ylabel('$\lambda I_\lambda(nW/m^2/sr)$','interpreter','latex','fontsize',18)
vline([wlmin wlmax],{'k','k'})

savename=strcat(savedir,'NISTspec');
print(savename,'-dpng');%close

if band==1
 meannist=[mean(nist_arr(1:2)), mean(nist_arr(3:4)),...
           mean(nist_arr(5:7)) mean(nist_arr(8:10))];
else
 meannist=[mean(nist_arr(1:3)), mean(nist_arr(4:6)),...
           mean(nist_arr(7:9)) mean(nist_arr(10:12))];
end
%% get imager ambient template 1
iset=1;
mapstack1=zeros(1024);maskstack1=zeros(1024);
for itime=1:numel(time_arr{iset})
time=char(time_arr{iset}(itime));
load(strcat(framedir,time,'_framedat'),'framedat');
map=framedat.map.*cal;
[~,mask]=get_skymap(map,mask_inst,3,3);maskmap=map.*mask;
mapstack1=mapstack1+maskmap;maskstack1=maskstack1+mask;
figure
imageclip(map);
title(strcat('TM',num2str(band),'\_ambient\_',time));
ylabel(colorbar, 'nW/m^2/sr','fontsize',15);
savename=strcat(savedir,'Imager',time);
print(savename,'-dpng');close
end
mapstack1=mapstack1./maskstack1;
mapstack1(find(mapstack1~=mapstack1))=0;
mapstack1(find(mapstack1==inf))=0;
mapstack1(find(mapstack1==-inf))=0;

%%% get imager ambient template 4
iset=4;
mapstack4=zeros(1024);maskstack4=zeros(1024);
for itime=1:numel(time_arr{iset})
time=char(time_arr{iset}(itime));
load(strcat(framedir,time,'_framedat'),'framedat');
map=framedat.map.*cal;
[~,mask]=get_skymap(map,mask_inst,3,3);maskmap=map.*mask;
mapstack4=mapstack4+maskmap;maskstack4=maskstack4+mask;
figure
imageclip(map);
title(strcat('TM',num2str(band),'\_ambient\_',time));
ylabel(colorbar, 'nW/m^2/sr','fontsize',15);
savename=strcat(savedir,'Imager',time);
print(savename,'-dpng');close
end
mapstack4=mapstack4./maskstack4;
mapstack4(find(mapstack4~=mapstack4))=0;
mapstack4(find(mapstack4==inf))=0;
mapstack4(find(mapstack4==-inf))=0;
%% get data sub ambient template
avg_calfac=[];std_calfac=[];
mapstack2=zeros(1024);maskstack2=zeros(1024);
mapstack3=zeros(1024);maskstack3=zeros(1024);
for iset=[2,3]
    if iset==2
        ambmap=mapstack1;ambmask=maskstack1;
    else
        ambmap=mapstack4;ambmask=maskstack4;
    end
for itime=1:numel(time_arr{iset})
time=char(time_arr{iset}(itime));
load(strcat(framedir,time,'_framedat'),'framedat');
map=framedat.map.*cal;
[~,mask]=get_skymap(map,mask_inst,3,3);maskmap=map.*mask;

if iset==2
mapstack2=mapstack2+maskmap;maskstack2=maskstack2+mask;
else
mapstack3=mapstack3+maskmap;maskstack3=maskstack3+mask;
end

ambmask(find(ambmask))=1;
submap=(map-ambmap).*mask.*ambmask;
[~,cmask]=get_skymap(submap,mask,3,3);submap=submap.*cmask;

submapff=submap./FF5;
submapff(find(submapff~=submapff))=0;
submapff(find(submapff==inf))=0;
submapff(find(submapff==-inf))=0;


figure
setwinsize(gcf,1000,1500)

subplot(2,3,1)
imageclip(map);
title('light map');
ylabel(colorbar, 'nW/m^2/sr','fontsize',15);

subplot(2,3,2)
imageclip(submap);
title('light map - ambient map');
ylabel(colorbar, 'nW/m^2/sr','fontsize',15);
v=caxis;

subplot(2,3,3)
imageclip(submapff);
title('(light map - ambient map)/FF');
caxis(v);

subplot(2,3,4)
h=histogram(submap(find(submap~=0)));hold on
histogram(submapff(find(submapff~=0)),'binedges',h.BinEdges);
legend({'lgiht-amb','(light-amb)/FF'})

calfac=submapff./meannist(iset);

calfac(find(calfac==inf))=0;calfac(find(calfac==-inf))=0;
subplot(2,3,5)
histogram(calfac(find(calfac~=0)),'Normalization','probability');
avgr=mean(calfac(find(calfac~=0)));
stdr=std(calfac(find(calfac~=0)));
title(sprintf('I_{CIBER}/I_{NIST}=%.2f +- %.4f',avgr,stdr),'fontsize',15);

figtitle(time,'fontweight','bold');

savename=strcat(savedir,'Imager',time);
print(savename,'-dpng');%close

avg_calfac=[avg_calfac avgr];
std_calfac=[std_calfac stdr];
end
end

save(strcat(savedir,'avg_calfac'),'avg_calfac');    
save(strcat(savedir,'std_calfac'),'std_calfac');
%%
figure;
errorbar(1:10,avg_calfac,std_calfac,'.k');
ylabel('$I_{\rm CIBER}/I_{\rm NIST}$','interpreter','latex','fontsize',18)
xlabel('exposure #','fontsize',18)
if band==1
    ylim([0.55 0.75]);
else
    ylim([0.35 0.65]);
end

xlim([0 11])
savename=strcat(savedir,'cal');
print(savename,'-dpng');%close
%%
maskjoint=ones(1024);
maskjoint(find(maskstack2==0))=0;maskjoint(find(maskstack3==0))=0;

mapstack2=mapstack2./maskstack2;
mapstack2(find(mapstack2~=mapstack2))=0;
mapstack2(find(mapstack2==inf))=0;mapstack2(find(mapstack2==-inf))=0;

mapstack3=mapstack3./maskstack3;
mapstack3(find(mapstack3~=mapstack3))=0;
mapstack3(find(mapstack3==inf))=0;mapstack3(find(mapstack3==-inf))=0;

ratiomap=(mapstack3./mapstack2).*maskjoint;
figure
setwinsize(gcf,1000,500)
subplot(2,2,1)
imageclip(ratiomap);
title('light sets ratio')

[~,ratiomask]=get_skymap(ratiomap,ones(1024),5,5);
ratiomap=ratiomap.*ratiomask;
subplot(2,2,2)
histogram(ratiomap(find(ratiomap~=0)));

maskjoint=ones(1024);
maskjoint(find(maskstack1==0))=0;maskjoint(find(maskstack4==0))=0;

mapstack1=mapstack1./maskstack2;
mapstack1(find(mapstack1~=mapstack1))=0;
mapstack1(find(mapstack1==inf))=0;mapstack1(find(mapstack1==-inf))=0;

mapstack4=mapstack4./maskstack4;
mapstack4(find(mapstack4~=mapstack4))=0;
mapstack4(find(mapstack4==inf))=0;mapstack4(find(mapstack4==-inf))=0;

ratiomap=(mapstack4./mapstack1).*maskjoint;
setwinsize(gcf,1000,500)
subplot(2,2,3)
imageclip(ratiomap);
title('ambient sets ratio')


[~,ratiomask]=get_skymap(ratiomap,ones(1024),5,5);
ratiomap=ratiomap.*ratiomask;
subplot(2,2,4)
histogram(ratiomap(find(ratiomap~=0)));
savename=strcat(savedir,'ratio');
print(savename,'-dpng');%close