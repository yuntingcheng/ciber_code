%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This code process the raw data to get FF and calmap, 
%   and save the processed data in alldat
% - The bigmask include Mike's mask, inst mask, 
%   crazy pixel mask, jack mask.
% - Unit: 
%   rawmap - ADU/fr
%   dsub,calmap: nW/m2/sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=36277;
band=2;
beam=7;
iter_mkk=10;
iter_clip=3;

cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

savedir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/';
wdir=strcat('/Users/ytcheng/ciber/doc/20161102_36277/NoiseModel/TM',...
                num2str(band),'/full/PS2D/');
%%% get inst_mask and DC template %%%%
DCdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/362xx/';
load(strcat(DCdir,'band',num2str(band),'_DCtemplate'),'DCtemplate');
load(strcat(DCdir,'band',num2str(band),'_mask_inst'),'mask_inst');
%%
%%% get field name and nfr
alldat(1).name = 'SWIRE';
alldat(2).name = 'NEP';
alldat(3).name = 'BootesB';
%%% load Mike's data
datadir=strcat('/Users/ytcheng/ciber/data/',num2str(flight),'/dr/');
rawdir=strcat('/Users/ytcheng/ciber/data/',num2str(flight),'/slopedata/');
for i=1:3
    dark_time=get_dark_times(flight,band,i);
    alldat(i).nfr=dark_time.nfr;
    
    dataname = strcat(datadir,'TM',num2str(band),'_',alldat(i).name,...
        '_dr130124.mat');load(dataname);
    load(sprintf('%sTM%d_%s',rawdir,band,alldat(i).name),'rawmap');
    alldat(i).rawmap=rawmap;%ADU/fr
    alldat(i).astrometry=data.astrometry;
    alldat(i).mzmask=double(~data.mask.mask);
    alldat(i).bl_l =data.psf.l;
    alldat(i).bl=data.psf.bl;
    alldat(i).dsub = (alldat(i).rawmap - DCtemplate).*cal; % nW/m2/sr
end
%% get bigmask: mzmask + crazy pixel mask + jack mask
%%% crazy pixel mask
for i=1:3
mzmask=alldat(i).mzmask;
[~,bigmask]=get_skymap(alldat(i).dsub,mzmask.*mask_inst,1);
alldat(i).bigmask=bigmask;
end

%%% jack mask
for i=1:3
load(sprintf('%sTM%d_%s_1st',rawdir,band,alldat(i).name),'rawmap1');
load(sprintf('%sTM%d_%s_2nd',rawdir,band,alldat(i).name),'rawmap2');
rawdiff=(rawmap1-rawmap2).*alldat(i).bigmask;
b=rawdiff(find(rawdiff~=0));
jackmask=ones(1024);
jackmask(find(rawdiff>mean(b)+2*std(b)))=0;
jackmask(find(rawdiff<mean(b)-2*std(b)))=0;
alldat(i).bigmask=jackmask.*alldat(i).bigmask;
end
%%
%%% stack FF
FF=zeros(1024);stack_mask=zeros(1024);
for i=1:3
bigmaski=alldat(i).bigmask;
obsi=alldat(i).dsub.*bigmaski;mean_obsi=mean(obsi(find(obsi)));
FF=FF+(obsi./sqrt(mean_obsi));
stack_mask=stack_mask+bigmaski.*sqrt(mean_obsi);
end
FF=FF./stack_mask;FF((find(FF~=FF)))=0;
save(strcat(savedir,'band',num2str(band),'_FF'),'FF');

%%% write bigmask include no flat pixels (make sure include crazy pix mask)
for i=1:3
bigmaski=alldat(i).bigmask;
bigmaski(find(stack_mask==0))=0;
alldat(i).bigmask=bigmaski;
end
%% get calmap
for i=1:3
dsub=alldat(i).dsub;
mask=alldat(i).bigmask;
calmap=dsub.*mask./FF;
calmap(find(calmap~=calmap))=0;
alldat(i).calmap=calmap;
end
save(strcat(savedir,'band',num2str(band),'_alldat'),'alldat');
%% get and save weighted Mkk & Mkkg for bigmask

load(strcat(savedir,'band',num2str(band),'_alldat'),'alldat');
for i=1:3
load(strcat(wdir,'b',num2str(band),'_i',...
      num2str(i), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';

disp('-------------------------------')
disp(i)
disp('-------------------------------')

[~,l,~,~,lbin]=get_angular_spec(alldat(i).dsub,alldat(i).dsub,beam);
maski=alldat(i).bigmask;
wMkk =get_mkk_sim(maski,beam,lbin,iter_mkk,numel(lbin),1,weight,0,NaN);
%wMkkg=get_mkk_grad(maski,beam,lbin,iter_mkk,1,weight);

alldat(i).wMkk=wMkk;
%alldat(i).wMkkg=wMkkg;
end
save(strcat(savedir,'band',num2str(band),'_alldat'),'alldat');