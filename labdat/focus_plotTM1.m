%%%%%%%%%%%%%plot TM1 nfr3,4,5%%%%%%%%%%%%%%%%%%
inst=1;
datadir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/fitTM'...
    ,num2str(inst),'/'); 
flight=40030;
cp=get_cal_params('flight',flight);
frate=cp(inst).framerate;%frame/sec
savedir='/Users/ytcheng/ciber/doc/20160912_CalFac/Focus';

figure;clf

nfr=3;
Tint=nfr/frate;
scanfile=dir(sprintf('%snfr%d*',datadir,nfr));

avg_arr=[];var_arr=[];
for i=1:numel(scanfile)
dataname=sprintf('%s%s',datadir,scanfile(i).name);
load(dataname,'fitdat');
avg=fitdat.avg_arr;var=fitdat.var_arr;
avg_arr=[avg_arr avg];var_arr=[var_arr var];  
end

pr=polyfit(var_arr,-avg_arr,1);
plt3=plot(-avg_arr,var_arr,'o','color',get_color(nfr));hold on
plot(polyval(pr,var_arr),var_arr,'color',get_color(nfr));
g1fit3=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);

nfr=4;
Tint=nfr/frate;
scanfile=dir(sprintf('%snfr%d*',datadir,nfr));

avg_arr=[];var_arr=[];
for i=1:numel(scanfile)
dataname=sprintf('%s%s',datadir,scanfile(i).name);
load(dataname,'fitdat');
avg=fitdat.avg_arr;var=fitdat.var_arr;
avg_arr=[avg_arr avg];var_arr=[var_arr var];  
end

pr=polyfit(var_arr,-avg_arr,1);
plt4=plot(-avg_arr,var_arr,'o','color',get_color(nfr));hold on
plot(polyval(pr,var_arr),var_arr,'color',get_color(nfr));
g1fit4=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);

nfr=5;
Tint=nfr/frate;
scanfile=dir(sprintf('%snfr%d*',datadir,nfr));

avg_arr=[];var_arr=[];
for i=1:numel(scanfile)
dataname=sprintf('%s%s',datadir,scanfile(i).name);
load(dataname,'fitdat');
avg=fitdat.avg_arr;var=fitdat.var_arr;
avg_arr=[avg_arr avg];var_arr=[var_arr var];  
end

pr=polyfit(var_arr,-avg_arr,1);
plt5=plot(-avg_arr,var_arr,'o','color',get_color(nfr));hold on
plot(polyval(pr,var_arr),var_arr,'color',get_color(nfr));
g1fit5=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);

xlabel('-mean (ADU/fr)','interpreter','latex','fontsize',18)
ylabel('$\sigma^2$(ADU/fr)','interpreter','latex','fontsize',18)

legend([plt3,plt4,plt5],{sprintf('3 frames, G1=%.3f',g1fit3),...
                         sprintf('4 frames, G1=%.3f',g1fit4),...
                         sprintf('5 frames, G1=%.3f',g1fit5)},...
    'Location','northwest','FontSize',15);
legend boxoff

imname=sprintf('%sTM%dall',savedir,inst);
print(imname,'-dpng');
%%
%%%%%%%%%%plot TM1 nfr 3,4,5 single plot%%%%%%%%%%%%%%
nfr=2;
inst=1;
datadir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/fitTM'...
    ,num2str(inst),'/'); 
flight=40030;
cp=get_cal_params('flight',flight);
frate=cp(inst).framerate;%frame/sec
savedir='/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/';

figure;clf

Tint=nfr/frate;
scanfile=dir(sprintf('%snfr%d*',datadir,nfr));

avg_arr=[];var_arr=[];
for i=1:numel(scanfile)
dataname=sprintf('%s%s',datadir,scanfile(i).name);
load(dataname,'fitdat');
avg=fitdat.avg_arr;var=fitdat.var_arr;
avg_arr=[avg_arr avg];var_arr=[var_arr var];  
end

pr=polyfit(var_arr,-avg_arr,1);
plt=plot(-avg_arr,var_arr,'o','color',get_color(1));hold on
plot(polyval(pr,var_arr),var_arr,'color',get_color(1));
g1fit=pr(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);


xlabel('-mean (ADU/fr)','interpreter','latex','fontsize',18)
ylabel('$\sigma^2$(ADU/fr)','interpreter','latex','fontsize',18)

legend([plt],{sprintf('%d frames, G1=%.3f',nfr,g1fit)},...
    'Location','southeast','FontSize',15);
legend boxoff

imname=sprintf('%sTM%dnfr%d',savedir,inst,nfr);
print(imname,'-dpng');
