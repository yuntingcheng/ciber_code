band=1;
flight=40030;
cp=get_cal_params('flight',flight);
frate=cp(band).framerate;%frame/sec

load(strcat('/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/band',...
    num2str(band),'_mask_inst'));
framedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/FF/TM'...
    ,num2str(band),'/framedat/'); 
savedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/FF/TM'...
    ,num2str(band),'/'); 

%%% get time labels%%%
%manually remove good frame<5 cases, see print out from G1_framefit.mat
switch band
case 1
%'10-33-57','10-41-26' removed
time1={{'10-17-17';'10-18-30';'10-19-07';'10-19-29'},...
    {'10-20-03';'10-20-51';'10-21-15'}};
time2={{'10-24-36';'10-25-40';'10-26-00';'10-26-21';'10-26-44'},...
    {'10-27-32';'10-27-51';'10-28-10';'10-28-50'}};
time3={{'10-33-00';'10-33-36';'10-34-18';'10-34-40'},...
    {'10-35-30';'10-35-53';'10-36-15'}};
time4={{'10-40-03';'10-40-24';'10-40-46'},...
    {'10-42-23';'10-42-49'}};
time_cell={time1,time2,time3,time4};
case 2
%'14-03-14' removed
time1={{'13-40-49';'13-41-48';'13-42-12';'13-42-31';...
        '13-43-08';'13-43-35';'13-44-04'},...
    {'13-44-40';'13-45-08';'13-45-30';'13-45-58';'13-46-19';'13-46-43'}};
time2={{'13-49-42';'13-50-18';'13-50-36';'13-50-54';...
        '13-51-13';'13-51-31';'13-52-08';'13-52-28';'13-52-47'},...
    {'13-53-37';'13-53-56';'13-54-15';'13-54-57';'13-55-32';'13-56-05';...
     '13-56-28';'13-57-03'}};
time3={{'13-58-24';'13-59-06';'13-59-42';'14-00-04';...
        '14-00-28';'14-00-51';'14-01-16';'14-01-41';'14-02-17'},...
    {'14-03-32';'14-03-51';'14-04-25';'14-04-44';'14-05-09';...
     '14-05-29';'14-06-13'}};
time_cell={time1,time2,time3};
end
%%
for nfr=[3,4,5]

%%%%%%%%%%% get ambient stack %%%%%%%%%%%%%
ambmap_arr=zeros(numel(time_cell),1024,1024);
ambmask_arr=zeros(numel(time_cell),1024,1024);

for iset=1:numel(time_cell)
mapstack=zeros(1024);maskstack=zeros(1024);    
for itime=1:numel(time_cell{iset}{2})
time=char(time_cell{iset}{2}{itime});
load(strcat(framedir,time,'_framedat'),'framedat');

map=eval(sprintf('framedat.map%d',nfr));maskframe=framedat.mask;
[~,mask]=get_skymap(map,maskframe,5,3);maskmap=map.*mask;

mapstack=mapstack+maskmap;maskstack=maskstack+mask;

figure
histogram(maskmap(find(maskmap~=0)));
title(time)
savename=strcat(savedir,'map/',time,'_hist',num2str(nfr));
%print(savename,'-dpng');close
figure
imageclip(maskmap);
title(time)
savename=strcat(savedir,'map/',time,'_map',num2str(nfr));
%print(savename,'-dpng');close

end
mapstack=mapstack./maskstack;
mapstack(find(mapstack~=mapstack))=0;
mapstack(find(mapstack==inf))=0;mapstack(find(mapstack==-inf))=0;
ambmap(iset,:,:)=mapstack;
maskstack(find(maskstack>0))=1;ambmask(iset,:,:)=maskstack;
end

%%%%%%%%%%% get light stack %%%%%%%%%%%%%
stackmap_arr=zeros(numel(time_cell),1024,1024);
stackmask_arr=zeros(numel(time_cell),1024,1024);

for iset=1:numel(time_cell)
amb=squeeze(ambmap(iset,:,:));mamb=squeeze(ambmask(iset,:,:));
mapstack=zeros(1024);maskstack=zeros(1024);    
for itime=1:numel(time_cell{iset}{1})
time=char(time_cell{iset}{1}{itime});
load(strcat(framedir,time,'_framedat'),'framedat');
map=eval(sprintf('framedat.map%d',nfr));maskframe=framedat.mask;

[~,mask]=get_skymap(map,maskframe.*mamb,5,3);maskmap=map.*mask;
mapstack=mapstack+maskmap;maskstack=maskstack+mask;

figure
histogram(maskmap(find(maskmap~=0)));
title(time)
savename=strcat(savedir,'map/',time,'_hist',num2str(nfr));
%print(savename,'-dpng');close
figure
imageclip(maskmap);
title(time)
savename=strcat(savedir,'map/',time,'_map',num2str(nfr));
%print(savename,'-dpng');close

end
mapstack=mapstack./maskstack;
mapstack(find(mapstack~=mapstack))=0;
mapstack(find(mapstack==inf))=0;mapstack(find(mapstack==-inf))=0;

maskstack(find(maskstack>0))=1;
[~,maskstack]=get_skymap(mapstack,maskstack,5,3);
mapstack=mapstack.*maskstack;

stackmap(iset,:,:)=mapstack;
stackmask(iset,:,:)=maskstack;
end
%%
%%% get diff variance 
var_arr=[];mean_arr=[];
for iset=1:numel(time_cell)
temp=squeeze(stackmap(iset,:,:));mtemp=squeeze(stackmask(iset,:,:));
for itime=1:numel(time_cell{iset}{1})-1
timei=char(time_cell{iset}{1}{itime});
load(strcat(framedir,timei,'_framedat'),'framedat');
map=eval(sprintf('framedat.map%d',nfr));maskframei=framedat.mask;
map=map.*maskframei;
[~,maski]=get_skymap(map,maskframe,5,3);maskmapi=map.*maski;

for jtime=itime+1:numel(time_cell{iset}{1})
timej=char(time_cell{iset}{1}{jtime});
load(strcat(framedir,timej,'_framedat'),'framedat');
map=eval(sprintf('framedat.map%d',nfr));maskframe=framedat.mask;
map=map.*maskframe;
[~,maskj]=get_skymap(map,maskframe,5,3);maskmapj=map.*maskj;

diffraw=(maskmapi-maskmapj)./sqrt(2).*maski.*maskj;

[~,maskd]=get_skymap(diffraw,maski.*maskj.*mask_inst,5,3);

%%% TM1 manual mask two strip 
if band==1;maskd(:,1:150)=0;maskd(:,500:650)=0;diff=diffraw.*maskd;end

%%% manual mask edge of bottom raw, left and right column
%%% do dc offset
if band==2;maskd(1:40,:)=0;maskd(:,1:40)=0;maskd(:,1000:end)=0;
diff=dc_offset_remove(diffraw,maskd);diff=diff.*maskd;end

%{
figure
imageclip(diffraw);
title(strcat(timei,'\_',timej))
savename=strcat(savedir,'diff/',timei,'_',timej,'_raw',num2str(nfr));
%print(savename,'-dpng');close

figure
histogram(diff(find(diff~=0)));
title(strcat(timei,'\_',timej))
savename=strcat(savedir,'diff/',timei,'_',timej,'_hist',num2str(nfr));
%print(savename,'-dpng');close
%}

figure
imageclip(diff);
title(strcat(timei,'\_',timej))
savename=strcat(savedir,'diff/',timei,'_',timej,'_map',num2str(nfr));
%print(savename,'-dpng');close


var_arr=[var_arr var(diff(find(diff~=0)))];
temp=temp.*maskd;
mean_arr=[mean_arr mean(temp(find(temp~=0)))];
end
end
end
%%
save(strcat(savedir,'var_arr',num2str(nfr)),'var_arr');
save(strcat(savedir,'mean_arr',num2str(nfr)),'mean_arr');
end
%%
figure
nfr=3;Tint=nfr/frate;c='b';
load(strcat(savedir,'var_arr',num2str(nfr)),'var_arr');
load(strcat(savedir,'mean_arr',num2str(nfr)),'mean_arr');
pr=polyfit(var_arr,-mean_arr,1);
g1fit3=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);
plt3=plot(-mean_arr,var_arr,strcat(c,'o'));hold on
plot(polyval(pr,var_arr),var_arr,c);

nfr=4;Tint=nfr/frate;c='r';
load(strcat(savedir,'var_arr',num2str(nfr)),'var_arr');
load(strcat(savedir,'mean_arr',num2str(nfr)),'mean_arr');
pr=polyfit(var_arr,-mean_arr,1);
g1fit4=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);
plt4=plot(-mean_arr,var_arr,strcat(c,'o'));hold on
plot(polyval(pr,var_arr),var_arr,c);

nfr=5;Tint=nfr/frate;c='k';
load(strcat(savedir,'var_arr',num2str(nfr)),'var_arr');
load(strcat(savedir,'mean_arr',num2str(nfr)),'mean_arr');
pr=polyfit(var_arr,-mean_arr,1);
g1fit5=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);
plt5=plot(-mean_arr,var_arr,strcat(c,'o'));hold on
plot(polyval(pr,var_arr),var_arr,c);

legend([plt3,plt4,plt5],{sprintf('3 frames, G1=%.3f',g1fit3),...
                         sprintf('4 frames, G1=%.3f',g1fit4),...
                         sprintf('5 frames, G1=%.3f',g1fit5)},...
    'Location','northwest','FontSize',15);
legend boxoff

xlabel('-mean (ADU/fr)','interpreter','latex','fontsize',18)
ylabel('$\sigma^2$(ADU/fr)','interpreter','latex','fontsize',18)

imname=strcat(savedir,'G1fit');
print(imname,'-dpng');

disp(sprintf('G1[3fr,4fr,5fr,mean]=[%.3f,%.3f,%.3f,%.3f]',...
            g1fit3,g1fit4,g1fit5,mean([g1fit3,g1fit4,g1fit5])));