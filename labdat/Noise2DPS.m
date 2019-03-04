%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This code calculate the read noise 2D PS full and diff with dark data.
% - The mask is the true sky mask (bigmask) from get_alldat.
% - Unit: nW/m2/sr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
band=1;
beam=7;
iter_mkk=10;
iter_clip=3; 

cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

savedir=strcat('/Users/ytcheng/ciber/doc/20160905_NoiseModel/TM',...
                    num2str(band),'/');
darkdir='/Users/ytcheng/ciber/doc/20150810_DarkCurrent/40030/';

DCdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(DCdir,'band',num2str(band),'_DCtemplate'),'DCtemplate');

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'));
load(strcat(alldatdir,'band',num2str(band),'_FF5'),'FF5');
%% get mean and std of nosie 2dPS of each field

for fieldnum=1:8
bigmask=alldat(fieldnum).bigmask;
fieldname=alldat(fieldnum).name;
nfr=alldat(fieldnum).nfr;

scandir=sprintf('%sTM%d/%s/full/data/',darkdir,band,fieldname);
scanfull=dir(strcat(scandir,'dark*'));
scandir=sprintf('%sTM%d/%s/first/data/',darkdir,band,fieldname);
scanfirst=dir(strcat(scandir,'dark*'));
scandir=sprintf('%sTM%d/%s/second/data/',darkdir,band,fieldname);
scansecond=dir(strcat(scandir,'dark*'));

PS2dfull_arr=zeros(numel(scanfull),1024,1024);
PS2ddiff_arr=zeros(numel(scanfull),1024,1024);
for i=1:numel(scanfull)
disp(sprintf('field=%d,count=%d/%d',fieldnum,i,numel(scanfull)));
%%%%%%%%%% full frames%%%%%%%%%%%%%%%%%%
fname=sprintf('%sTM%d/%s/full/data/%s',...
        darkdir,band,fieldname,scanfull(i).name);
dark=load(fname);dark=dark.darkmap;
dark=dark-DCtemplate;dark=dark.*cal;dark=dark./FF5;
dark(find(dark~=dark))=0;dark(find(dark==inf|dark==-inf))=0;
[~,mask]=get_skymap(dark,bigmask,iter_clip);

if band==2;dark=dc_offset_remove(dark,mask);end

dark=dark.*mask;meandark=mean(dark(find(dark~=0)));
dark=(dark-meandark).*mask;
[~,~,~,~,lbin,~,nCl2d] = get_angular_spec(dark,dark,beam);
PS2dfull_arr(i,:,:)=nCl2d;
%%%%%%%%%% diff %%%%%%%%%%%%%%%%%%%%%%%%
fname1=sprintf('%sTM%d/%s/first/data/%s',...
        darkdir,band,fieldname,scanfirst(i).name);
dark1=load(fname1);dark1=dark1.darkmap;
fname2=sprintf('%sTM%d/%s/second/data/%s',...
        darkdir,band,fieldname,scansecond(i).name);
dark2=load(fname2);dark2=dark2.darkmap;
diff=(dark1-dark2)./sqrt(2);diff=diff.*cal;diff=diff./FF5;
diff(find(diff~=diff))=0;diff(find(diff==inf|diff==-inf))=0;
[~,mask]=get_skymap(diff,bigmask,iter_clip);

if band==2;diff=dc_offset_remove(diff,mask);end

diff=diff.*mask;meandiff=mean(diff(find(diff~=0)));
diff=(diff-meandiff).*mask;
[~,~,~,~,lbin,~,nCl2d] = get_angular_spec(diff,diff,beam);
PS2ddiff_arr(i,:,:)=nCl2d;
imageclip(nCl2d);
end

save(strcat(savedir,'full/PS2D/b',num2str(band),'_i', ...
      num2str(fieldnum), '_PS2d_arr'),'PS2dfull_arr');
save(strcat(savedir,'diff/PS2D/b',num2str(band),'_i', ...
      num2str(fieldnum), '_PS2d_arr'),'PS2ddiff_arr');
  
std_noise=squeeze(std(PS2dfull_arr));
save(strcat(savedir,'full/PS2D/b',num2str(band),'_i', ...
      num2str(fieldnum), '_std_noise'),'std_noise');

mean_noise=squeeze(mean(PS2dfull_arr));
save(strcat(savedir,'full/PS2D/b',num2str(band),'_i', ...
      num2str(fieldnum), '_mean_noise'),'mean_noise');
end

