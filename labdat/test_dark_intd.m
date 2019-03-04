%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Test sigma_dark map - N frames relation
% - use the diff to eliminate dark current
% - Ref: Garnett & Forrest 1993 Eq 19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band=1;

darkdir=strcat('/Users/ytcheng/ciber/data/DarkDat/06-05-2013/TM',...
    num2str(band),'/');
savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_dark_int/';
iter_clip=3;

switch band
case 1
time_arr={'13-53-13';'13-53-13';'13-54-36';'14-23-09';'14-23-09';...
          '14-41-29';'14-41-29';'14-45-38';'14-45-38';'14-49-42';...
          '14-49-42';'15-01-10';'15-01-10';'15-10-18';'15-10-18';...
          '15-26-56';'15-29-12';'15-29-12'};
start_arr=[3,21,3,3,55,3,49,19,68,3,28,3,33,3,51,13,3,66];
end_arr=[18,34,30,52,79,46,78,65,130,22,83,30,71,48,102,71,63,105];
case 2
time_arr={'13-53-13';'13-53-13';'13-54-36';'14-23-09';'14-23-09';...
          '14-41-29';'14-41-29';'14-45-38';'14-45-38';'14-49-42';...
          '15-01-10';'15-01-10';'15-10-18';'15-10-18';'15-26-56';...
          '15-29-12';'15-29-12'};
start_arr=[3,21,3,3,58,3,50,19,68,3,3,59,13,78,3,3,28];
end_arr=[18,33,30,55,79,47,78,65,131,23,56,70,75,102,36,25,90];
end

dur_arr=end_arr-start_arr+1;
[~,sortind]=sort(dur_arr);
% remove the 9th data in band1
if band==1;sortind=sortind(sortind~=9);end
%% linefit dark map fr by fr
%{
for iset=1:floor(numel(sortind)/2)
iset1=sortind(iset);
iset2=sortind(iset+1);
nfr_arr=2:3:dur_arr(iset1);
sig_arr=zeros(1,numel(nfr_arr));
diffmap_arr=zeros(numel(nfr_arr),1024,1024);
diffmask_arr=zeros(numel(nfr_arr),1024,1024);

time1=time_arr{iset1};    
scanfile1=dir(strcat(darkdir,'*',time1,'*'));
time2=time_arr{iset2};    
scanfile2=dir(strcat(darkdir,'*',time2,'*'));

for i=1:numel(nfr_arr)
nfr=nfr_arr(i);
startfr1=start_arr(iset1);
endfr1=startfr1+nfr-1;
startfr2=start_arr(iset2);
endfr2=startfr2+nfr-1;

frames1=zeros(nfr,1024,1024);
for j=startfr1:endfr1
    fname=strcat(darkdir,scanfile1(j).name);
    frame = imrotate(fitsread(fname),270);
    frames1(j-startfr1+1,:,:)=frame;
end
darkmap1=fastlinefit_frin(frames1,0);

frames2=zeros(nfr,1024,1024);
for j=startfr2:endfr2
    fname=strcat(darkdir,scanfile2(j).name);
    frame = imrotate(fitsread(fname),270);
    frames2(j-startfr2+1,:,:)=frame;
end
darkmap2=fastlinefit_frin(frames2,0);

diffmap=(darkmap1-darkmap2)./sqrt(2);

mask=mask_clip(diffmap,iter_clip);
sig=std(diffmap(find(mask)));
sig_arr(i)=sig;
diffmap_arr(i,:,:)=diffmap;
diffmask_arr(i,:,:)=mask;

imageclip(diffmap);

pr=sprintf('%d/%d,nfr=%d,time1=%s,time2=%s,sig=%.2e',...
               i,numel(nfr_arr),nfr,time1,time2,sig);disp(pr);
end

darktest(iset).time1=time1;
darktest(iset).startfr1=start_arr(iset1);
darktest(iset).endfr1=end_arr(iset1);
darktest(iset).time2=time2;
darktest(iset).startfr2=start_arr(iset2);
darktest(iset).endfr2=end_arr(iset2);
darktest(iset).nfr_arr=nfr_arr;
darktest(iset).sig_arr=sig_arr;
darktest(iset).diffmap_arr=diffmap_arr;
darktest(iset).diffmask_arr=diffmask_arr;
end
save(strcat(savedir,'band',num2str(band),'_darktestdiff'),'darktest');
%}
%% fit the sigma_RN and sigma_DC
load(strcat(savedir,'band',num2str(band),'_darktestdiff'),'darktest');

figure
for i=1:numel(darktest)
    
nfr_arr=darktest(i).nfr_arr;
sig_arr=darktest(i).sig_arr;

if numel(nfr_arr)>numel(nfrplt_arr);nfrplt_arr=nfr_arr;end

plot(nfr_arr,sig_arr,'-o','color',[0.6 0.6 0.6]);hold on
end
xlabel('$N$','interpreter','latex','fontsize',18)
ylabel('$\sigma[ADU/fr]$','interpreter','latex','fontsize',18)

cp=get_cal_params('flight',40030);
frate=cp(band).framerate;%frame/sec

if band==1;sig_cds=10;elseif band==2;sig_cds=9;end
savename=strcat(savedir,'band',num2str(band),'_var_N_diff_org');
print(savename,'-dpng');



figure
ydata_arr=[];
xdata_arr=[];
nfrplt_arr=[];
for i=1:numel(darktest)
    
nfr_arr=darktest(i).nfr_arr;
sig_arr=darktest(i).sig_arr;

if numel(nfr_arr)>numel(nfrplt_arr);nfrplt_arr=nfr_arr;end

xdata_arr=[xdata_arr nfr_arr];
ydata_arr=[ydata_arr sqrt(sig_arr.^2.*nfr_arr.*(nfr_arr.^2-1)./6)];


plot(nfr_arr,sqrt(sig_arr.^2.*nfr_arr.*(nfr_arr.^2-1)./6),...
        '-o','color',[0.6 0.6 0.6]);hold on
end
xlabel('$N$','interpreter','latex','fontsize',18)
ylabel('$\sqrt{\sigma^2N(N^2-1)/6}[ADU/fr]$',...
    'interpreter','latex','fontsize',18)

cp=get_cal_params('flight',40030);
frate=cp(band).framerate;%frame/sec

if band==1;sig_cds=10;elseif band==2;sig_cds=9;end
pr=sprintf('TM%d, sig_CDS/dt=%.3f+-%.3f[ADU/fr],sig_CDS/dt=%.3f[e-/s]',...
            band,mean(ydata_arr),std(ydata_arr),sig_cds*frate);disp(pr);
pr=sprintf('TM%d,G1 fit=%.3f+-%.3f, G1=%.3f',...
     band,sig_cds*frate/mean(ydata_arr),...
sig_cds*frate*std(ydata_arr)/mean(ydata_arr)^2,cp(band).apf2eps);disp(pr);

savename=strcat(savedir,'band',num2str(band),'_var_N_diff');
print(savename,'-dpng');
