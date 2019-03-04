
%%%%%%%%%%plot TM2 nfr 3,4,5 single plot%%%%%%%%%%%%%%
nfr=5;
inst=2;
datadir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/fitTM'...
    ,num2str(inst),'/'); 
flight=40030;
cp=get_cal_params('flight',flight);
frate=cp(inst).framerate;%frame/sec
savedir='/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/';

figure;clf

Tint=nfr/frate;
scanfile=dir(sprintf('%snfr%d*',datadir,nfr));

set1=3:10;set2=1:17;set2(set1)=[];

avg1_arr=[];var1_arr=[];
for i=set1
dataname=sprintf('%s%s',datadir,scanfile(i).name);
load(dataname,'fitdat');
avg=fitdat.avg_arr;var=fitdat.var_arr;
avg1_arr=[avg1_arr avg];var1_arr=[var1_arr var]; 
end
pr1=polyfit(var1_arr,-avg1_arr,1);
plt1=plot(-avg1_arr,var1_arr,'o','color',get_color(1));hold on
plot(polyval(pr1,var1_arr),var1_arr,'color',get_color(1));
g1fit1=pr1(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);

avg2_arr=[];var2_arr=[];
for i=set2
dataname=sprintf('%s%s',datadir,scanfile(i).name);
load(dataname,'fitdat');
avg=fitdat.avg_arr;var=fitdat.var_arr;
avg2_arr=[avg2_arr avg];var2_arr=[var2_arr var]; 
end
pr2=polyfit(var2_arr,-avg2_arr,1);
plt2=plot(-avg2_arr,var2_arr,'o','color',get_color(2));hold on
plot(polyval(pr2,var2_arr),var2_arr,'color',get_color(2));
g1fit2=pr2(1)*6*(nfr^2+1)/Tint/5/(nfr^2-1);

xlabel('-mean (ADU/fr)','interpreter','latex','fontsize',18)
ylabel('$\sigma^2$(ADU/fr)','interpreter','latex','fontsize',18)

legend([plt1,plt2],{sprintf('%d frames, G1=%.3f',nfr,g1fit1),...
                    sprintf('%d frames, G1=%.3f',nfr,g1fit2) },...
    'Location','southeast','FontSize',15);
legend boxoff

imname=sprintf('%sTM%dnfr%d',savedir,inst,nfr);
print(imname,'-dpng');
