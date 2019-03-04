savedir = '/Users/ytcheng/Desktop/NEP/';

band=1;
flight=36277;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

alldatdir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'));

map = alldat(2).rawmap.*cal;
calmap = alldat(2).calmap;
ra = alldat(2).astrometry.ra;
dec = alldat(2).astrometry.dec;
mask = alldat(2).bigmask;


fits_write(strcat(savedir,'NEP36277_map_I.fits'),map);
fits_write(strcat(savedir,'NEP36277_calmap_I.fits'),calmap);
fits_write(strcat(savedir,'NEP36277_mask_I.fits'),mask);
fits_write(strcat(savedir,'NEP36277_ra_I.fits'),ra);
fits_write(strcat(savedir,'NEP36277_dec_I.fits'),dec);
%%
flight=36265;
band=2;
iter_clip=3;

cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

load(strcat('/Volumes/HD1TB/CIBER/data/',num2str(flight),'/dr/TM',...
    num2str(band),'_NEP_dr130124.mat'));

map = data.release.all.map;
mask = double(~data.release.all.mask);
ra = data.release.all.astrometry.ra;
dec = data.release.all.astrometry.dec;

fits_write(strcat(savedir,'NEP36265_map_H.fits'),map);
fits_write(strcat(savedir,'NEP36265_mask_H.fits'),mask);
fits_write(strcat(savedir,'NEP36265_ra_H.fits'),ra);
fits_write(strcat(savedir,'NEP36265_dec_H.fits'),dec);
%%
savedir = '/Users/ytcheng/Desktop/NEP/';

band=2;
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'));

map = alldat(2).rawmap.*cal;
calmap = alldat(2).calmap5;
ra = alldat(2).astrometry.ra;
dec = alldat(2).astrometry.dec;
mask = alldat(2).bigmask;


fits_write(strcat(savedir,'NEP40030_map_H.fits'),map);
fits_write(strcat(savedir,'NEP40030_calmap_H.fits'),calmap);
fits_write(strcat(savedir,'NEP40030_mask_H.fits'),mask);
fits_write(strcat(savedir,'NEP40030_ra_H.fits'),ra);
fits_write(strcat(savedir,'NEP40030_dec_H.fits'),dec);
%% 40030 DGL
savedir = '/Users/ytcheng/Desktop/DGL/';

band=2; name = 'H';
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'));

map = alldat(1).rawmap.*cal;
calmap = alldat(1).calmap5;
ra = alldat(1).astrometry.ra;
dec = alldat(1).astrometry.dec;
mask = alldat(1).bigmask;


fits_write(strcat(savedir,'DGL40030_map_',name,'.fits'),map);
fits_write(strcat(savedir,'DGL40030_calmap_',name,'.fits'),calmap);
fits_write(strcat(savedir,'DGL40030_mask_',name,'.fits'),mask);
fits_write(strcat(savedir,'DGL40030_ra_',name,'.fits'),ra);
fits_write(strcat(savedir,'DGL40030_dec_',name,'.fits'),dec);
