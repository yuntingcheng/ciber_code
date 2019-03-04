%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given flight and band, go through the full and diff dark and 
%take the PS before and after Mkk. Save the PS in struc darkps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
    
pixsize=7;
iter_clip=3;
%%
darkdir='/Users/ytcheng/ciber/doc/20150810_DarkCurrent/';
savedir=strcat('/Users/ytcheng/ciber/doc/20160906_NoiseRealization/',...
            'darkstat/',num2str(flight),'/');
alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
        
fields=get_fields(flight,inst);

%%
for ifield=1:numel(fields)
fieldname=fields(ifield).name;

fulldir=sprintf('%s/%d/TM%d/%s/full/data/',darkdir,flight,inst,fieldname);
scanfull=dir(strcat(fulldir,'dark*'));

ndarks=floor(numel(scanfull)/2);
Clf_arr=zeros(ndarks,29);
Cl2d_arr=zeros(ndarks,1024,1024);

bigmask=alldat(ifield).bigmask;
mClf_arr=zeros(ndarks,29);
mCl2d_arr=zeros(ndarks,1024,1024);

for i=1:2:2*ndarks
load(strcat(fulldir,scanfull(i).name));dark1=darkmap;
load(strcat(fulldir,scanfull(i+1).name));dark2=darkmap;

dark=(dark1-dark2)./sqrt(2);
[~,mdark]=get_skymap(dark,ones(1024),iter_clip);
dark=dc_offset_remove(dark,mdark);
dark=dark-mean(dark(find(dark)));dark=dark.*mdark;
[Cl,l,~,~,~,~,Cl2d]=get_angular_spec(dark,dark,pixsize);
Clf_arr((i+1)/2,:)=Cl;
Cl2d_arr((i+1)/2,:,:)=Cl2d;

dark=(dark1-dark2)./sqrt(2);
[~,mdark]=get_skymap(dark,bigmask,iter_clip);
dark=dc_offset_remove(dark,mdark);
dark=dark-mean(dark(find(dark)));dark=dark.*mdark;
[mCl,l,~,~,~,~,mCl2d]=get_angular_spec(dark,dark,pixsize);
mClf_arr((i+1)/2,:)=mCl;
mCl2d_arr((i+1)/2,:,:)=mCl2d;

disp(sprintf('%d,TM%d,%s,%d/%d full done',...
                flight,inst,fieldname,i,numel(scanfull)));

end
darkps(ifield).flight=flight;
darkps(ifield).inst=inst;
darkps(ifield).field=fieldname;
darkps(ifield).l=l;
darkps(ifield).Clf_arr=Clf_arr;
darkps(ifield).Clf2d_arr=Cl2d_arr;
darkps(ifield).Clf2d_ave=squeeze(mean(Cl2d_arr));
darkps(ifield).Clf2d_std=squeeze(std(Cl2d_arr));

mdarkps(ifield).flight=flight;
mdarkps(ifield).inst=inst;
mdarkps(ifield).field=fieldname;
mdarkps(ifield).l=l;
mdarkps(ifield).mClf_arr=mClf_arr;
mdarkps(ifield).mClf2d_arr=mCl2d_arr;
mdarkps(ifield).mClf2d_ave=squeeze(mean(mCl2d_arr));
mdarkps(ifield).mClf2d_std=squeeze(std(mCl2d_arr));
end
%%
%%%%%%%%%%%%% half %%%%%%%%%%%%%%%%
for ifield=1:numel(fields)
fieldname=fields(ifield).name;

dir1=sprintf('%s/%d/TM%d/%s/first/data/',darkdir,flight,inst,fieldname);
dir2=sprintf('%s/%d/TM%d/%s/second/data/',darkdir,flight,inst,fieldname);

scan1=dir(strcat(dir1,'dark*'));
scan2=dir(strcat(dir2,'dark*'));

Cld_arr=zeros(numel(scan1),29);
Cl2d_arr=zeros(numel(scan1),1024,1024);

bigmask=alldat(ifield).bigmask;
mCld_arr=zeros(numel(scan1),29);
mCl2d_arr=zeros(numel(scan1),1024,1024);

for i=1:numel(scan1)
load(strcat(dir1,scan1(i).name));dark1=darkmap;
load(strcat(dir2,scan2(i).name));dark2=darkmap;

dark=(dark1-dark2)./sqrt(2);
[~,mdark]=get_skymap(dark,ones(1024),iter_clip);
dark=dc_offset_remove(dark,mdark);
dark=dark-mean(dark(find(dark)));dark=dark.*mdark;
[Cl,l,~,~,~,~,Cl2d]=get_angular_spec(dark,dark,pixsize);
Cld_arr(i,:)=Cl;
Cl2d_arr(i,:,:)=Cl2d;

dark=(dark1-dark2)./sqrt(2);
[~,mdark]=get_skymap(dark,bigmask,iter_clip);
dark=dc_offset_remove(dark,mdark);
dark=dark-mean(dark(find(dark)));dark=dark.*mdark;
[mCl,l,~,~,~,~,mCl2d]=get_angular_spec(dark,dark,pixsize);
mCld_arr(i,:)=mCl;
mCl2d_arr(i,:,:)=mCl2d;

disp(sprintf('%d,TM%d,%s,%d/%d diff done',...
                flight,inst,fieldname,i,numel(scanfull)));

end
darkps(ifield).Cld_arr=Cld_arr;
darkps(ifield).Cld2d_arr=Cl2d_arr;
darkps(ifield).Cld2d_ave=squeeze(mean(Cl2d_arr));
darkps(ifield).Cld2d_std=squeeze(std(Cl2d_arr));
mdarkps(ifield).mCld_arr=mCld_arr;
mdarkps(ifield).mCld2d_arr=mCl2d_arr;
mdarkps(ifield).mCld2d_ave=squeeze(mean(mCl2d_arr));
mdarkps(ifield).mCld2d_std=squeeze(std(mCl2d_arr));

end
%%
save(strcat(savedir,'TM',num2str(inst),'_darkps'),'darkps');
save(strcat(savedir,'TM',num2str(inst),'_mdarkps'),'mdarkps');
