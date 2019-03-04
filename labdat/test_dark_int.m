%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Test sigma_dark map - N frames relation
% - Ref: Garnett & Forrest 1993 Eq 19
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band=2;

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
%% linefit dark map fr by fr
%{
for iset=1:numel(time_arr)
time=time_arr{iset};    
scanfile=dir(strcat(darkdir,'*',time,'*'));
nfr_arr=2:3:dur_arr(iset);
sig_arr=zeros(1,numel(nfr_arr));

for i=1:numel(nfr_arr)
nfr=nfr_arr(i);
startfr=start_arr(iset);
endfr=startfr+nfr-1;

frames=zeros(nfr,1024,1024);
for j=startfr:endfr
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    frames(j-startfr+1,:,:)=frame;
end
darkmap=fastlinefit_frin(frames,0);

mask=mask_clip(darkmap,iter_clip);
sig=std(darkmap(find(mask)));
sig_arr(i)=sig;

pr=sprintf('%d/%d,%s,nfr=%d,sig=%.2e',i,numel(nfr_arr),time,nfr,sig);
disp(pr);
end

darktest(iset).time=time;
darktest(iset).startfr=start_arr(iset);
darktest(iset).endfr=end_arr(iset);
darktest(iset).dur=dur_arr(iset);
darktest(iset).nfr_arr=nfr_arr;
darktest(iset).sig_arr=sig_arr;
darktest(iset).darkmap=darkmap;%the longest integration time darkmap
end
save(strcat(savedir,'band',num2str(band),'_darktest'),'darktest');
%}
%% fit the sigma_RN and sigma_DC
load(strcat(savedir,'band',num2str(band),'_darktest'),'darktest');

figure
ydata_arr=[];
xdata_arr=[];
nfrplt_arr=[];
for i=1:numel(darktest)
    
nfr_arr=darktest(i).nfr_arr;
sig_arr=darktest(i).sig_arr;

if numel(nfr_arr)>numel(nfrplt_arr);nfrplt_arr=nfr_arr;end

if ~(band==1 & i==9)
% TM1 9th(14-45-38 68-130) sig^2 goes very high at large N, 
% don't use in fitting
xdata_arr=[xdata_arr nfr_arr];
ydata_arr=[ydata_arr sig_arr.^2.*nfr_arr.*(nfr_arr.^2-1)];
end

plot(nfr_arr,sig_arr.^2.*nfr_arr.*(nfr_arr.^2-1),...
        '-o','color',[0.6 0.6 0.6]);hold on
end

% do the parameter fit
% initial value y(1)=180,y(2)=2e-3 is arbitrary chosen by eye fitting
F = @(y,x)y(1)+y(2).*x.*(x.^2-1);
[fitpars] = lsqcurvefit(F,[200 1e-3],xdata_arr,ydata_arr);

fit=plot(nfrplt_arr,fitpars(1)+...
    fitpars(2).*nfrplt_arr.*(nfrplt_arr.^2-1),'k','linewidth',2);


cp=get_cal_params('flight',40030);
frate=cp(band).framerate;%frame/sec

if band==1;sig_cds=10;elseif band==2;sig_cds=9;end
pr=sprintf('TM%d, sig_CDS/dt=%.3f[ADU/fr],sig_CDS/dt=%.3f[e-/s]',...
            band,sqrt(fitpars(1)/6),sig_cds*frate);disp(pr);
pr=sprintf('TM%d,G1 fit=%.3f, G1=%.3f',...
     band,sig_cds*frate/sqrt(fitpars(1)/6),cp(band).apf2eps);disp(pr);

xlabel('$N$','interpreter','latex','fontsize',18)
ylabel('$\sigma^2N(N^2-1)[(ADU/fr)^2]$',...
    'interpreter','latex','fontsize',18)

legend([fit],{sprintf('y=%.2f+%.4f N(N^2-1)',...
    fitpars(1),fitpars(2))},'Location','northwest','FontSize',15);
legend boxoff
savename=strcat(savedir,'band',num2str(band),'_var_N');
print(savename,'-dpng');
