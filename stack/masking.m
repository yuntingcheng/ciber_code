flight=40030;
mypaths=get_paths(flight);
inst=1;
pixscale=7;

cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

quad_arr=['A','B','C','D'];
srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
load(strcat(mypaths.alldat,'TM',num2str(inst),'/','maskdat'),'maskdat');
load(strcat(mypaths.alldat,'TM',num2str(inst),'/','FFdat'),'FFdat');

load(sprintf('%sTM%d/darklongdat',mypaths.filtmap,inst),'darklongdat');
DCtemplate=darklongdat.DCtemplate; clear darklongdat
%%
ifield=4;
dt=get_dark_times(flight,inst,ifield);

FF=FFdat(ifield).FF;
%%% get flight map %%%
loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
load(strcat(loaddir,'flightmap'),'flightmap');
rawmap=flightmap.filtmapf;
calmap=(rawmap-DCtemplate)./FF;
calmap(find(calmap~=calmap))=0;
calmap(find(calmap==-inf))=0;
calmap(find(calmap==inf))=0;

%%% get sim srcmap %%%
ukmap=zeros(1024);
tmmap=zeros(1024);
for iquad=1:4
    quad=quad_arr(iquad);
    sukmap=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt.fits'));
    stmmap=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapt_2m.fits'));
    
    if iquad==1
        ukmap(1:512,1:512)=sukmap;
        tmmap(1:512,1:512)=stmmap;
    elseif iquad==2
        ukmap(513:1024,1:512)=sukmap;
        tmmap(513:1024,1:512)=stmmap;
    elseif iquad==3
        ukmap(1:512,513:1024)=sukmap;
        tmmap(1:512,513:1024)=stmmap;
    else
        ukmap(513:1024,513:1024)=sukmap;
        tmmap(513:1024,513:1024)=stmmap;
    end
end
%%% get masks %%%%
bigmask=maskdat.mask(ifield).bigmask;
nosrcmask=maskdat.mask(ifield).nosrc;
%% plot to see the map
figure
setwinsize(gcf,1500,500)
subplot(1,3,1)
imageclip(log10(abs(calmap.*cal)));
title('CIBER')
subplot(1,3,2)
imageclip(log10(abs(ukmap)));
title('UKIDSS sim map')
subplot(1,3,3)
imageclip(log10(abs(tmmap)));
title('2MASS sim map')
%% get bigmask Cl
[~,~,~,~,binl]=get_angular_spec(randn(1024),randn(1024),pixscale);

[~,mask1]=get_skymap(calmap,bigmask,5,5);
%[cCl,wCl,Cl2d,l,binl,dCl]=get_Cl(map,mask,Mkk,beam,weight,varargin);
mkk=numel(mask1)./numel(find(mask1));

cbmap=calmap.*cal./rcal;
cbmap1=(cbmap-mean(cbmap(find(mask1)))).*mask1;
ukmap1=(ukmap-mean(ukmap(find(mask1)))).*mask1;
[cCl,l]=get_angular_spec(cbmap1,cbmap1,pixscale);
[uCl,l]=get_angular_spec(ukmap1,ukmap1,pixscale);
[xCl,l]=get_angular_spec(cbmap1,ukmap1,pixscale);
cCl=cCl.*mkk;
uCl=uCl.*mkk;
xCl=xCl.*mkk;
%% Masking test: change alpha
alpha_arr=[-10,-15,-20];
figure
setwinsize(gcf,1000,1000)
for i=1:numel(alpha_arr)
m_max=17;
alpha=alpha_arr(i);
beta=-m_max*alpha+3;

strmask1=make_strmask(flight,inst,ifield,'2m',alpha,beta,m_max);
mask1=nosrcmask.*strmask1;
[~,mask1]=get_skymap(calmap,mask1,5,5);
[~,mask1]=get_skymap(ukmap,mask1,5,5);

mkk=numel(mask1)./numel(find(mask1));

subplot(2,2,i)
imageclip(mask1.*calmap);
title(strcat('\alpha=',num2str(alpha)));

cbmap=calmap.*cal;
cbmap1=(cbmap-mean(cbmap(find(mask1)))).*mask1;
ukmap1=(ukmap-mean(ukmap(find(mask1)))).*mask1;
[cCl,l]=get_angular_spec(cbmap1,cbmap1,pixscale);
[uCl,l]=get_angular_spec(ukmap1,ukmap1,pixscale);
[xCl,l]=get_angular_spec(cbmap1,ukmap1,pixscale);
cCl=cCl.*mkk;
uCl=uCl.*mkk;
xCl=xCl.*mkk;

subplot(2,2,4)
loglog(l,l.*(l+1).*cCl./2./pi,'-','color',get_color(i),...
    'Displayname',strcat('CBxCB, \alpha=',num2str(alpha)));hold on
loglog(l,l.*(l+1).*uCl./2./pi,'--','color',get_color(i),...
    'Displayname',strcat('UKxUK, \alpha=',num2str(alpha)));
loglog(l,l.*(l+1).*xCl./2./pi,':','color',get_color(i),...
    'Displayname',strcat('CBxUK, \alpha=',num2str(alpha)));
end
xlim([1e2,2e5]);
title(dt.name);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend('show','location','northwest')
legend boxoff

%% Masking test: change depth
figure
m_max_arr=[17,18,19];
setwinsize(gcf,1000,1000)
for i=1:numel(alpha_arr)
m_max=m_max_arr(i);
alpha=-15;
beta=-17*alpha+3;

strmask1=make_strmask(flight,inst,ifield,'2m',alpha,beta,m_max);
strmask2=make_strmask(flight,inst,ifield,'uk',alpha,beta,m_max);

mask1=nosrcmask.*strmask1.*strmask2;
[~,mask1]=get_skymap(calmap,mask1,5,5);
[~,mask1]=get_skymap(ukmap,mask1,5,5);
mkk=numel(mask1)./numel(find(mask1));

subplot(2,2,i)
imageclip(mask1.*calmap);
title(strcat('mAB=',num2str(m_max)));

cbmap=calmap.*cal;
cbmap1=(cbmap-mean(cbmap(find(mask1)))).*mask1;
ukmap1=(ukmap-mean(ukmap(find(mask1)))).*mask1;
[cCl,l]=get_angular_spec(cbmap1,cbmap1,pixscale);
[uCl,l]=get_angular_spec(ukmap1,ukmap1,pixscale);
[xCl,l]=get_angular_spec(cbmap1,ukmap1,pixscale);
cCl=cCl.*mkk;
uCl=uCl.*mkk;
xCl=xCl.*mkk;

subplot(2,2,4)
loglog(l,l.*(l+1).*cCl./2./pi,'-','color',get_color(i),...
    'Displayname',strcat('CBxCB, m_{AB}=',num2str(m_max)));hold on
loglog(l,l.*(l+1).*uCl./2./pi,'--','color',get_color(i),...
    'Displayname',strcat('UKxUK, m_{AB}=',num2str(m_max)));
loglog(l,l.*(l+1).*xCl./2./pi,':','color',get_color(i),...
    'Displayname',strcat('CBxUK, m_{AB}=',num2str(m_max)));
end
xlim([1e2,2e5]);
title(dt.name);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend('show','location','northwest')
legend boxoff