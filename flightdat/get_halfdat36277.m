%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This code gives the 1st and 2nd raw and calmap, and their diff.
% - The bigmask,wMkk,wMkkg same as alldat
% - Unit: 
%   rawfirst/rawsecond - ADU/fr
%   others maps and means: nW/m2/sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=36277;
band=2;
beam=7;
iter_clip=3;

cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

rawdir=strcat('/Users/ytcheng/ciber/data/',num2str(flight),'/slopedata/');
alldatdir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'),'alldat');
load(strcat(alldatdir,'band',num2str(band),'_FF'),'FF');
DCdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/362xx/';
load(strcat(DCdir,'band',num2str(band),'_DCtemplate'),'DCtemplate');

savedir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDiff/';
%%
%%% read 1st and 2nd halves raw data
for i=1:3
load(sprintf('%sTM%d_%s_1st',rawdir,band,alldat(i).name),'rawmap1');
load(sprintf('%sTM%d_%s_2nd',rawdir,band,alldat(i).name),'rawmap2');
halfdat(i).name=alldat(i).name;
halfdat(i).rawfirst=rawmap1;
halfdat(i).rawsecond=rawmap2;
halfdat(i).nfr=floor(alldat(i).nfr/2);
halfdat(i).bigmask=alldat(i).bigmask;
halfdat(i).wMkk=alldat(i).wMkk;
%halfdat(i).wMkkg=alldat(i).wMkkg;

end
clear alldat
%%
for i=1:3
bigmask=halfdat(i).bigmask;
raw1=(halfdat(i).rawfirst).*cal;raw2=(halfdat(i).rawsecond).*cal;
rawd=(raw1-raw2)./sqrt(2);

cal1=(raw1-(DCtemplate.*cal)).*bigmask./FF;cal1(find(cal1~=cal1))=0;
cal2=(raw2-(DCtemplate.*cal)).*bigmask./FF;cal2(find(cal2~=cal2))=0;
cald=(cal1-cal2)./sqrt(2);

halfdat(i).raw1=raw1.*bigmask;halfdat(i).raw2=raw2.*bigmask;
halfdat(i).rawd=rawd.*bigmask;

halfdat(i).cal1=cal1.*bigmask;halfdat(i).cal2=cal2.*bigmask;
halfdat(i).cald=cald.*bigmask;


halfdat(i).mean1=mean(cal1(find(cal1)));
halfdat(i).mean2=mean(cal2(find(cal2)));

figure;
imageclip(cal1);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_cal1');
title(halfdat(i).name);print(imname,'-dpng');close

figure;
imageclip(cal2);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_cal2');
title(halfdat(i).name);print(imname,'-dpng');close

figure;
imageclip(cald);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_cald');
title(halfdat(i).name);print(imname,'-dpng');close

end
save(strcat(savedir,'band',num2str(band),'_halfdat'),'halfdat');