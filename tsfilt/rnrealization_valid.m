function rnrealization_valid(flight,inst,ifield)

pixscale=7;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);


loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(savedir,'darkstat'),'darkstat');

load(sprintf('%sTM%d/darklongdat',mypaths.filtmap,inst),'darklongdat');
DCtemplate=darklongdat.DCtemplate_nomask; clear darklongdat

load(strcat(savedir,'maskdat'),'maskdat');
load(strcat(savedir,'fwdat'),'fwdat');
bigmask=maskdat.mask(ifield).bigmask;
fw=fwdat(ifield).fw_filt;


%% get lab noise data

Cl_arr=zeros(numel(dt.time),29);
wCl_arr=zeros(numel(dt.time),29);

disp(sprintf('get lab RN PS'));

for i=1:numel(dt.time)
    load(strcat(loaddir,'labmap',num2str(i)),'labmap');
    filtmapf=labmap.filtmapf-DCtemplate;
    [~,maskf]=get_skymap(filtmapf,bigmask,4,5);
    filtmapf=filtmapf-mean(filtmapf(find(maskf)));filtmapf=filtmapf.*maskf;
    filtmapf=dc_offset_remove(filtmapf,maskf).*maskf;
    [Cl,l]=get_angular_spec(filtmapf,filtmapf,pixscale);
    [wCl]=get_angular_spec(filtmapf,filtmapf,pixscale,'w',fw);
    
    Cl_arr(i,:)=Cl;
    wCl_arr(i,:)=wCl;
end
Cl_arr=Cl_arr(:,9:29);
wCl_arr=wCl_arr(:,9:29);
%% noise realization

disp(sprintf('get sim RN PS'));

nsim=50;
simCl_arr=zeros(nsim,29);
simwCl_arr=zeros(nsim,29);
for i=1:nsim
    rnmap=readnoise_realization(darkstat(ifield).Cl2d_std,pixscale,'norand',1);
    
    [~,mask]=get_skymap(rnmap,bigmask,4,5);
    rnmap=rnmap-mean(rnmap(find(mask)));rnmap=rnmap.*mask;

    [Cl]=get_angular_spec(rnmap,rnmap,pixscale);
    [wCl]=get_angular_spec(rnmap,rnmap,pixscale,'w',fw);
    
    simCl_arr(i,:)=Cl;
    simwCl_arr(i,:)=wCl;
end
simCl_arr=simCl_arr(:,9:29);
simwCl_arr=simwCl_arr(:,9:29);
%% knox for comparison
[cCl,l,~,dCl]=get_angular_spec(randn(1024),randn(1024),pixscale);
l=l(find(cCl));cCl=cCl(find(cCl));dCl=dCl(find(dCl));
knoxratio=dCl./cCl;
%% make plots
fig=figure;
setwinsize(gcf,1000,500)

subplot(1,2,1)
yvalue=(l.*(l+1).*(mean(Cl_arr))./2./pi);
ebar=(l.*(l+1).*(mean(Cl_arr).*knoxratio)./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltknox=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',10);hold on

y1=(l.*(l+1).*(prctile(Cl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(Cl_arr,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

y1=(l.*(l+1).*(prctile(simCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(simCl_arr,84))./2./pi);
pltsim=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [1,0.6,0.8],'facealpha',0.5,'EdgeColor','none');hold on


xlim([1e2,2e5]);ylim([1e-6,1e-2]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title(strcat(dt.name,'--no Fweight'));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi$','interpreter','latex','fontsize',18)
legend([pltknox,pltn,pltsim],...
 {'data-knox','dark data','RN sim'},'Location','southeast','FontSize',15);
legend boxoff

subplot(1,2,2)
yvalue=(l.*(l+1).*(mean(wCl_arr))./2./pi);
ebar=(l.*(l+1).*(mean(wCl_arr).*knoxratio)./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltknox=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',10);hold on

y1=(l.*(l+1).*(prctile(wCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(wCl_arr,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

y1=(l.*(l+1).*(prctile(simwCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(simwCl_arr,84))./2./pi);
pltsim=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [1,0.6,0.8],'facealpha',0.5,'EdgeColor','none');hold on


xlim([1e2,2e5]);ylim([1e-6,1e-2]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title(strcat(dt.name,'--Fweight'));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi$','interpreter','latex','fontsize',18)
legend([pltknox,pltn,pltsim],...
 {'data-knox','dark data','RN sim'},'Location','southeast','FontSize',15);
legend boxoff

return
