%%
flight=40030;
inst=2;
mypaths=get_paths(flight);
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

savedir='/Users/ytcheng/ciber/doc/20170904_External/plots/';
% get the dark current
loaddir=sprintf('%sTM%d/',mypaths.filtmap,inst);
load(sprintf('%s/darklongdat',loaddir),'darklongdat');
DCtemplate=darklongdat.DCtemplate; clear darklongdat

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/FFdat',loaddir),'FFdat');
load(sprintf('%s/maskdat',loaddir),'maskdat');

zldir=strcat('/Users/ytcheng/ciber/doc/20170904_External/ZL/kelsall/TM',...
    num2str(inst),'/');

% ZL scaling from Kesall 1.25 to CB band
%%%% to be updated !!!!
if inst==1
    rzl=1.165;
    rcal=0.58;
else
    rzl=0.538;
    rcal=0.38;
end

% get psf fit dat
psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

%%
zlmean_arr=[];
cbmz_arr=[];
cbyt_arr=[];

for ifield=4:8
    disp(sprintf('ifield=%d',ifield));
    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'flightmap'),'flightmap');
    dt=get_dark_times(flight,inst,ifield);

    bigmask=maskdat.mask(ifield).bigmask;
    FF=FFdat(ifield).FF;
    
    flightf=flightmap.rawmapf;
    flightf=(flightf-DCtemplate)./FF.*bigmask;
    flightf(find(flightf~=flightf))=0;
    flightf(find(flightf==-inf))=0;
    flightf(find(flightf==inf))=0;
    cbmean=mean(flightf(find(flightf)));
    cbmz_arr=[cbmz_arr cbmean.*cal];
    
    corr=fitpsfdat(ifield).bestparam_norm(1)+fitpsfdat(ifield).bestparam_norm(2);
    corr=corr*100;
    cbyt_arr=[cbyt_arr cbmean.*cal./rcal.*corr];
    
    % get Kelsall
    zlkmap=fitsread(strcat(zldir,dt.name,'_band1.fits'));
    zlkmap=zlkmap.*rzl;
    zlmean=mean(zlkmap(:));
    zlmean_arr=[zlmean_arr zlmean];
end
%%
plot(zlmean_arr,cbmz_arr,'bo-');hold on
plot(zlmean_arr,cbyt_arr,'ro-');
plot([1,2000],[1,2000],'k--');
if inst==1
    ylim([250,1100])
    xlim([250,800])
else
    xlim([100,400])
    ylim([100,400])
end
title(strcat('TM',num2str(inst)),'fontsize',15)
xlabel('Kelsall mean $nW/m^2/sr$','interpreter','latex','fontsize',18)
ylabel('CIBER mean $nW/m^2/sr$','interpreter','latex','fontsize',18)
legend({'MZ cal','YT cal','Kelsall=CIBER'},'fontsize',15,'location','northwest')
legend boxoff
imname=strcat(savedir,'TM',num2str(inst),'_Kelsall_CB_mean');
%print(imname,'-dpng');%close