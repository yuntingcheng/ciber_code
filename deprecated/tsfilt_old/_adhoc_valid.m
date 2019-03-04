%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - make the darkstat with only sigclip mask. 
% - adhoc noise validation for full and diff.
% - Conclusion: is not very consistent, especially for diff case.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1=-3.6;
load('/Volumes/HD1TB/CIBER/tsfilt/sinfiltamp/DarkLong/DCtemplate',...
                                                    'DCtemplate');
FFdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(FFdir,'band',num2str(inst),'_FF5'),'FF5');

savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
savedirhalf='/Users/ytcheng/ciber/doc/20170209_TsFilter/halfPS/';
load(strcat(savedirfull,'fwfulldat'),'fwfulldat');
load(strcat(savedirhalf,'fwhalfdat'),'fwhalfdat');
savedir='/Users/ytcheng/ciber/doc/20170320_FFsim/';
%% get and save dark stat
%{
for ifield=[8,7,6,5,4,2,1]

disp(sprintf('ifield=%d',ifield));
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
load(strcat(loaddir,'maskin'),'maskin');
dt=get_dark_times(flight,inst,ifield);

disp(sprintf('stack Cl2d_ave half'));
Cl2dd_ave=zeros(1024);
for i=1:numel(dt.time)
    load(strcat(loaddir,'labdat',num2str(i)),'labdat');
    
    filtmapd=(labdat.filtmap1-labdat.filtmap2)./2;
    [~,maskd]=get_skymap(filtmapd,ones(1024),4,5);
    [~,~,~,~,~,~,Cl2d]=get_angular_spec...
                    (filtmapd.*maskd,filtmapd.*maskd,pixscale);
    Cl2dd_ave=Cl2dd_ave+Cl2d./numel(dt.time);
end

disp(sprintf('stack Cl2d_ave full'));
Cl2df_ave=zeros(1024);
for i=2:2:numel(dt.time)
    load(strcat(loaddir,'labdat',num2str(i-1)),'labdat');
    map1=labdat.filtmapf;
    load(strcat(loaddir,'labdat',num2str(i)),'labdat');
    map2=labdat.filtmapf;
    filtmapd=(map1-map2)./sqrt(2);
    [~,maskd]=get_skymap(filtmapd,ones(1024),4,5);
    [~,~,~,~,~,~,Cl2d]=get_angular_spec...
                    (filtmapd.*maskd,filtmapd.*maskd,pixscale);
    Cl2df_ave=Cl2df_ave+Cl2d./numel(2:2:numel(dt.time));
end

darkstat(ifield).Cl2dd_ave=Cl2dd_ave;
darkstat(ifield).Cl2df_ave=Cl2df_ave;

end
save(strcat(savedir,'darkstat'),'darkstat');
%}
%%
load(strcat(savedir,'darkstat'),'darkstat');
%%% get Knox
[cCl,l,~,dCl]=get_angular_spec(randn(1024),randn(1024),pixscale);
l=l(find(cCl));cCl=cCl(find(cCl));dCl=dCl(find(dCl));
knoxratio=dCl./cCl;

ifield=1;
disp(sprintf('ifield=%d',ifield));
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
load(strcat(loaddir,'maskin'),'maskin');
dt=get_dark_times(flight,inst,ifield);
%% get lab noise data

disp(sprintf('get weighted Cl'));

nClf_arr=zeros(numel(dt.time),21);
nCld_arr=zeros(numel(dt.time),21);
for i=1:numel(dt.time)
    load(strcat(loaddir,'labdat',num2str(i)),'labdat');
    
    filtmapf=labdat.filtmapf;
    filtmapf=(filtmapf-DCtemplate); 
    W=(fftshift(fftshift(1./fwfulldat(ifield).std_fCl2d)))';
    [Clf]=get_Cl(filtmapf,maskin,fwfulldat(ifield).Mkk,pixscale,W);
    nClf_arr(i,:)=Clf;
    
    filtmapd=(labdat.filtmap1-labdat.filtmap2)./2;
    [~,maskd]=get_skymap(filtmapd,maskin,4,5);
    W=(fftshift(fftshift(1./fwhalfdat(ifield).std_dCl2d)))';
    [Cld]=get_Cl(filtmapd,maskd,fwhalfdat(ifield).Mkk,pixscale,W);
    
    nCld_arr(i,:)=Cld;
end
%% noise realization
nsim=50;
Wf=(fftshift(fftshift(1./fwfulldat(ifield).std_fCl2d)))';
Wh=(fftshift(fftshift(1./fwhalfdat(ifield).std_dCl2d)))';
Mkkf=fwfulldat(ifield).Mkk;
Mkkh=fwhalfdat(ifield).Mkk;
simClf_arr=zeros(nsim,21);
simCld_arr=zeros(nsim,21);
for i=1:nsim
    rnmap=readnoise_realization(darkstat(ifield).Cl2df_ave,...
                        pixscale,'norand',1);
    [Clf]=get_Cl(rnmap,maskin,Mkkf,pixscale,Wf);
    simClf_arr(i,:)=Clf;
    
    rnmap=readnoise_realization(darkstat(ifield).Cl2dd_ave,...
                        pixscale,'norand',1);
    [Cld]=get_Cl(rnmap,maskin,Mkkh,pixscale,Wh);
    simCld_arr(i,:)=Cld;
end
%%
fig=figure;
setwinsize(gcf,1000,500)

%%% plot PS full
subplot(1,2,1)
yvalue=(l.*(l+1).*(mean(nClf_arr)).*cal.^2./2./pi);
ebar=(l.*(l+1).*(mean(nClf_arr).*knoxratio)*cal.^2./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltknox=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',10);hold on

y1=(l.*(l+1).*(prctile(nClf_arr.*cal.^2,16))./2./pi);
y2=(l.*(l+1).*(prctile(nClf_arr.*cal.^2,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

y1=(l.*(l+1).*(prctile(simClf_arr.*cal.^2,16))./2./pi);
y2=(l.*(l+1).*(prctile(simClf_arr.*cal.^2,84))./2./pi);
pltsim=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [1,0.6,0.8],'facealpha',0.5,'EdgeColor','none');hold on


xlim([1e2,2e5]);ylim([1e-2,1e4]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title(strcat(dt.name,'--full'));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltknox,pltn,pltsim],...
 {'data-knox','dark data','RN sim'},...
       'Location','southeast','FontSize',15);
legend boxoff

%%% plot PS diff
subplot(1,2,2)
yvalue=(l.*(l+1).*(mean(nCld_arr)).*cal.^2./2./pi);
ebar=(l.*(l+1).*(mean(nCld_arr).*knoxratio)*cal.^2./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltknox=errorbar(l,yvalue,ebarlow,ebar,'.k','markersize',10);hold on

y1=(l.*(l+1).*(prctile(nCld_arr.*cal.^2,16))./2./pi);
y2=(l.*(l+1).*(prctile(nCld_arr.*cal.^2,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

y1=(l.*(l+1).*(prctile(simCld_arr.*cal.^2,16))./2./pi);
y2=(l.*(l+1).*(prctile(simCld_arr.*cal.^2,84))./2./pi);
pltsim=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [1,0.6,0.8],'facealpha',0.5,'EdgeColor','none');hold on


xlim([1e2,2e5]);ylim([1e-2,1e4]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title(strcat(dt.name,'--half'));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltknox,pltn,pltsim],...
 {'data-knox','dark data','RN sim'},...
       'Location','southeast','FontSize',15);
legend boxoff
imname=strcat(savedirhalf,'b',num2str(inst),'_i',num2str(ifield));
%print(imname,'-dpng');close



