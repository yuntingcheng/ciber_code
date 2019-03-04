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

savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/cats/';
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
%% Masking test: catalogs
alpha=-15;
m_max=17;
beta=250;

pscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','PSC');
xscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSC');
xscrmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSCrej');
pscrmask1=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','A');
pscrmask2=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','B');
pscrmask3=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','C');
pscrmask4=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','D');
pscrmask5=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','E');
pscrmask6=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
    'catname','PSCrej','rel','F');

strmask=pscmask.*xscmask.*xscmask.*pscrmask1;
%% clean the mask
maskin=strmask.*nosrcmask;
[~,maskin]=get_skymap(flightmap.rawmapf,maskin,4,5);

[~,maskin]=get_skymap(flightmap.filtmapf,maskin,4,5);

%%%%%%% jack mask %%%%%%
 diffmap=flightmap.rawmap1-flightmap.rawmap2;
[~,jackmask]=get_skymap(diffmap,maskin,3,5);

%%%% some hand mask determined by eye(ex:Cosmic Rays) %%%%%%
if inst==1 & ifield==8
jackmask=circular_mask(442,40,40,jackmask);
jackmask=circular_mask(742,342,10,jackmask);
end

if inst==2 & ifield==8
jackmask=circular_mask(450,395,10,jackmask);
jackmask=circular_mask(0,800,60,jackmask);
end

if inst==2 & ifield==5
jackmask=elliptical_mask(215,220,60,15,80,jackmask);
jackmask=elliptical_mask(710,90,80,20,60,jackmask);
end


[~,mask,meanmap,~,~,clipmax,clipmin]=get_skymap(calmap.*cal,jackmask,4,5);
mask(find(ukmap.*mask > meanmap+clipmax))=0;
mask(find(tmmap.*mask > meanmap+clipmax))=0;
%%
figure
setwinsize(gcf,1000,1000)

subplot(2,2,1)
imageclip(log10(calmap.*cal.*mask));
title('CIBER elat10')
[sm,fillmap]=fillpadsmooth(calmap.*cal,mask,10);
subplot(2,2,2)
imageclip(sm);
title('CIBER elat10 smooth')


subplot(2,2,3)
imageclip(log10(ukmap.*mask));
title('UKsim elat10')
[sm,fillmap]=fillpadsmooth(ukmap,mask,10);
subplot(2,2,4)
imageclip(sm);
title('UKsim elat10 smooth')

savename=strcat(savedir,'smooth');
print(savename,'-dpng');