inst=1;
pixscale=7;
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1=-3.6;G2=cal./G1;

savedir='/Users/ytcheng/ciber/doc/20160921_Simulation/';
weightdir=strcat('/Users/ytcheng/ciber/doc/20160905_NoiseModel/TM',...
                    num2str(inst),'/full/PS2D/');
alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');

drdir=strcat('/Users/ytcheng/ciber/data/',num2str(flight),'/dr/');
dataname = strcat(drdir,'TM',num2str(inst),'_',alldat(1).name,...
    '_dr150206.mat');% lab flat are the same for all fields
load(dataname);labflat=data.flat.labflat;

zldir = '/Users/ytcheng/ciber/data/TM1_zl_pmk/';
%% simulation input info
% input signal, extropolate Zemcov 2013 Fig S13
[~,l,~,~,lbin]=get_angular_spec(randn(1024),randn(1024),pixscale);
logx=[log10(2e3) log10(1e2)];x=10.^logx;
logy=[log10(2) log10(12)];y=10.^logy;
logCl=interp1(logx,logy,log10(l),'linear','extrap');
Clshot=ones(size(logCl)).*1e3.*2.*pi./1e5./(1e5+1);
Clin=(((10.^logCl)).*2.*pi./l./(l+1))+Clshot;
for i=1:8
    load(strcat(weightdir,'b',num2str(inst),'_i',...
      num2str(i), '_std_noise'),'std_noise');
    weight=(fftshift(fftshift(1./std_noise)))';

    name=alldat(i).name;
    nfr=alldat(i).nfr;

    zlname = strcat(zldir,'TM1_40030_zlgrad',alldat(i).name,'.fits');
    zlgrad = fliplr(fitsread_orient(zlname));

    meanmap=alldat(i).calmap5;
    meanmap=mean(meanmap(find(meanmap)));
    
    sigmap = map_from_power_spec(lbin,Clin,1024,1024,pixscale,1);

    siminfo(i).name=name;
    siminfo(i).nfr=nfr;
    siminfo(i).weight=weight;
    siminfo(i).meanmap=meanmap;
    siminfo(i).bigmask=alldat(i).bigmask;
    siminfo(i).zlgrad=zlgrad;
    siminfo(i).Mkk=alldat(i).wMkk;
    siminfo(i).sigmap=sigmap;
    %siminfo(i).Mkkg=alldat(i).wMkkg;
end
clear alldat;
%% 100 runs
tic
goods=[0,0,0,1,1,1,1,1];
for runs=1:100
    
%%%%% get sim maps %%%%%%
for i =1:8
[~,rnmap,~]=noise_realization(flight,inst,i,0,G1);
slopemap=(sigmap+siminfo(i).zlgrad).*labflat./cal;%% wrong!!! siminfo(i).sigmap
[phmap]=photonnoise_realization(slopemap,G1,siminfo(i).nfr,frate);
mapdat(i).sigmap=siminfo(i).sigmap;
mapdat(i).phmap=phmap.*cal;
mapdat(i).rnmap=rnmap.*cal;
end

%%%%% stack FF %%%%%
FF=zeros(size(labflat));stack_mask=zeros(size(labflat));
    
for i=1:8
bigmask=siminfo(i).bigmask;
zlgrad=siminfo(i).zlgrad;
sigmap=mapdat(i).sigmap;
phmap=mapdat(i).phmap;
rnmap=mapdat(i).rnmap;


obs=((sigmap+zlgrad).*labflat)+phmap+rnmap;
obs=obs-mean(obs(:))+siminfo(i).meanmap;
obs=obs.*bigmask;
mapdat(i).obsmap=obs;

if goods(i)
mean_obs=mean(obs(find(obs)));
FF=FF+(obs./sqrt(mean_obs));
stack_mask=stack_mask+bigmask.*sqrt(mean_obs);
end
end
FF=FF./stack_mask;FF((find(FF~=FF)))=0;

%%%%% get PS and save %%%%%
for i=1:8
bigmask=siminfo(i).bigmask;
zlgrad=siminfo(i).zlgrad;
weight=siminfo(i).weight;
Mkk=siminfo(i).Mkk;

sigmap=mapdat(i).sigmap;
phmap=mapdat(i).phmap;
rnmap=mapdat(i).rnmap;
frnmap=rnmap./FF;frnmap(find(frnmap~=frnmap))=0;

[cCl0]=get_Cl(sigmap+phmap+frnmap,bigmask,Mkk,pixscale,weight);
[cCls]=get_Cl(sigmap,bigmask,Mkk,pixscale,weight);


obs=mapdat(i).obsmap;
map1=obs./FF;map1(find(map1~=map1))=0;map1=map1.*bigmask;
%plane=plane_fit(map1,bigmask);
%map2=(map1-plane).*bigmask;

[cCl1]=get_Cl(map1,bigmask,Mkk,pixscale,weight);
%[cCl2]=get_Cl(map2,bigmask,Mkkg,beam,weight);

[cCln]=get_Cl(frnmap+phmap,bigmask,Mkk,pixscale,weight);

Cldat(i).cCl0=cCl0;Cldat(i).cCls=cCls;
Cldat(i).Cl1=cCl1;%Cldat(i).Cl2=cCl2;
Cldat(i).Cln=cCln;    
end

save(strcat(savedir,'Cldatruns/run',num2str(runs),'_Cldat'),'Cldat');
pr=sprintf('run=%d,t=%.1f min',runs,toc/60);disp(pr);
end
%%
%%%%%%%%% get PS data arr%%%%%%%%%%%%%%%%%%
for i=1:8
Cl0_arr=zeros(100,21);Cls_arr=zeros(100,21);
Cl1_arr=zeros(100,21);rCl1_arr=zeros(100,21);
%Cl2_arr=zeros(100,29);rCl2_arr=zeros(100,29);
Cln_arr=zeros(100,21);rCln_arr=zeros(100,21);
subCl0_arr=zeros(100,21);
subCl1_arr=zeros(100,21);

bigmask=siminfo(i).bigmask;
sigmap=mapdat(i).sigmap;
Mkk=siminfo(i).Mkk;
weight=siminfo(i).weight;
[cCls,~,~,l]=get_Cl(sigmap,bigmask,Mkk,pixscale,weight);

for runs=1:100
load(strcat(savedir,'Cldatruns/run',num2str(runs),'_Cldat'),'Cldat');
Cl0_arr(runs,:)=Cldat(i).cCl0;

Cls_arr(runs,:)=cCls;
Cl1_arr(runs,:)=Cldat(i).Cl1;
%Cl2_arr(runs,:)=Cldat(i).Cl2;
Cln_arr(runs,:)=Cldat(i).Cln;
subCl0_arr(runs,:)=Cldat(i).cCl0-Cldat(i).Cln;
subCl1_arr(runs,:)=cCls-Cldat(i).Cln;

rCl1_arr(runs,:)=Cldat(i).Cl1./Cldat(i).cCl0;
%rCl2_arr(runs,:)=Cldat(i).Cl2./Cldat(i).Cl0;
rCln_arr(runs,:)=Cldat(i).Cln./Cldat(i).cCl0;
end
save(sprintf('%sCldat/TM%d/i%d_Cl0_arr',savedir,inst,i),'Cl0_arr');
save(sprintf('%sCldat/TM%d/i%d_Cls_arr',savedir,inst,i),'Cls_arr');
save(sprintf('%sCldat/TM%d/i%d_Cl1_arr',savedir,inst,i),'Cl1_arr');
%save(sprintf('%sCldat/TM%d/i%d_Cl2_arr',savedir,inst,i),'Cl2_arr');
save(sprintf('%sCldat/TM%d/i%d_Cln_arr',savedir,inst,i),'Cln_arr');
save(sprintf('%sCldat/TM%d/i%d_rCl1_arr',savedir,inst,i),'rCl1_arr');
%save(sprintf('%sCldat/TM%d/i%d_rCl2_arr',savedir,inst,i),'rCl2_arr');
save(sprintf('%sCldat/TM%d/i%d_rCln_arr',savedir,inst,i),'rCln_arr');
save(sprintf('%sCldat/TM%d/i%d_subCl0_arr',savedir,inst,i),'subCl0_arr');
save(sprintf('%sCldat/TM%d/i%d_subCl1_arr',savedir,inst,i),'subCl1_arr');
end

% now clear the data in Cldatruns
%if ~exist(strcat(savedir,'Cldatruns/'), 'dir');mkdir(foldname);end
%%
for i=1:8
load(sprintf('%sCldat/TM%d/i%d_Cl0_arr',savedir,inst,i),'Cl0_arr');
load(sprintf('%sCldat/TM%d/i%d_Cls_arr',savedir,inst,i),'Cls_arr');
load(sprintf('%sCldat/TM%d/i%d_Cl1_arr',savedir,inst,i),'Cl1_arr');
%load(sprintf('%sCldat/TM%d/i%d_Cl2_arr',savedir,inst,i),'Cl2_arr');
load(sprintf('%sCldat/TM%d/i%d_Cln_arr',savedir,inst,i),'Cln_arr');
load(sprintf('%sCldat/TM%d/i%d_rCl1_arr',savedir,inst,i),'rCl1_arr');
%load(sprintf('%sCldat/TM%d/i%d_rCl2_arr',savedir,inst,i),'rCl2_arr');
load(sprintf('%sCldat/TM%d/i%d_rCln_arr',savedir,inst,i),'rCln_arr');
load(sprintf('%sCldat/TM%d/i%d_subCl0_arr',savedir,inst,i),'subCl0_arr');
load(sprintf('%sCldat/TM%d/i%d_subCl1_arr',savedir,inst,i),'subCl1_arr');

fig=figure;
pltCls=plot(l,l.*(l+1).*prctile(Cls_arr,50)./2./pi,...
                    '.','color',[0,0,0],'markersize',10);hold on
pltCl0=errorbar(l,l.*(l+1).*prctile(Cl0_arr,50)./2./pi,...
     l.*(l+1).*(prctile(Cl0_arr,50)-prctile(Cl0_arr,16))./2./pi,...
     l.*(l+1).*(prctile(Cl0_arr,84)-prctile(Cl0_arr,50))./2./pi,...
                    '.','color',[0.5,0.5,0.5],'markersize',10);hold on
pltCl1=errorbar(l,l.*(l+1).*prctile(Cl1_arr,50)./2./pi,...
     l.*(l+1).*(prctile(Cl1_arr,50)-prctile(Cl1_arr,16))./2./pi,...
     l.*(l+1).*(prctile(Cl1_arr,84)-prctile(Cl1_arr,50))./2./pi,...
                    '.r','markersize',10);hold on

y1=(l.*(l+1).*(prctile(Cln_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(Cln_arr,84))./2./pi);
pltCln=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on


xlim([2e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltCls,pltCl0,pltCl1,pltCln],...
 {'input signal','signal+noise map','output','noise model'},...
       'Location','southeast','FontSize',15);
legend boxoff
title(siminfo(i).name); 
savename=strcat(savedir,'plots/b',num2str(inst),'i',num2str(i),'_simPS');
print(savename,'-dpng');%close

end
%%
for i=1:8
load(sprintf('%sCldat/TM%d/i%d_Cl0_arr',savedir,inst,i),'Cl0_arr');
load(sprintf('%sCldat/TM%d/i%d_Cls_arr',savedir,inst,i),'Cls_arr');
load(sprintf('%sCldat/TM%d/i%d_Cl1_arr',savedir,inst,i),'Cl1_arr');
%load(sprintf('%sCldat/TM%d/i%d_Cl2_arr',savedir,inst,i),'Cl2_arr');
load(sprintf('%sCldat/TM%d/i%d_Cln_arr',savedir,inst,i),'Cln_arr');
load(sprintf('%sCldat/TM%d/i%d_rCl1_arr',savedir,inst,i),'rCl1_arr');
%load(sprintf('%sCldat/TM%d/i%d_rCl2_arr',savedir,inst,i),'rCl2_arr');
load(sprintf('%sCldat/TM%d/i%d_rCln_arr',savedir,inst,i),'rCln_arr');
load(sprintf('%sCldat/TM%d/i%d_subCl0_arr',savedir,inst,i),'subCl0_arr');
load(sprintf('%sCldat/TM%d/i%d_subCl1_arr',savedir,inst,i),'subCl1_arr');

fig=figure;
pltCls=plot(l,l.*(l+1).*prctile(Cls_arr,50)./2./pi,...
                    '.','color',[0,0,0],'markersize',20);hold on               
pltCl0=errorbar(l,l.*(l+1).*(mean(Cl0_arr)-mean(Cln_arr))./2./pi,...
     l.*(l+1).*(sqrt(std(Cl0_arr).^2+std(Cln_arr).^2))./2./pi,...
                    '.','color',[0.5,0.5,0.5],'markersize',10);hold on
pltCl1=errorbar(l,l.*(l+1).*(mean(Cl1_arr)-mean(Cln_arr))./2./pi,...
     l.*(l+1).*(sqrt(std(Cl1_arr).^2+std(Cln_arr).^2))./2./pi,...
                    'r.','markersize',10);hold on

y1=(l.*(l+1).*(prctile(Cln_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(Cln_arr,84))./2./pi);
pltCln=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on


xlim([1e2,2e5]);ylim([1e-2,5e3]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltCls,pltCl0,pltCl1,pltCln],...
 {'input signal','signal+noise map - noise model','output-noise model','noise model'},...
       'Location','southeast','FontSize',15);
legend boxoff
title(siminfo(i).name); 
savename=strcat(savedir,'plots/b',num2str(inst),'i',num2str(i),'_simPSsub');
print(savename,'-dpng');%close

end
