flight=40030;
inst=1;
ifield=8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

%% masking parameters
alpha=-9.7*1.2;
beta=159*1.2;
m_arr=10:1:21;

for m_max=m_arr
    
pscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','PSC');
xscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSC');
xscrmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSCrej');
pscrmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
'catname','PSCrej','rel','A');
ukmask=make_strmask_uk(flight,inst,ifield,alpha,beta,m_max);
strmask=pscmask.*xscmask.*xscrmask.*pscrmask.*ukmask;

fits_write(strcat(savedir,'masks/',dt.name,'_strmask',num2str(m_max)),strmask);
disp(sprintf('strmask within m=%d save',m_max));
end
%% get inst mask and map
mask=fits_read(strcat(savedir,'masks/',dt.name,'_mask.fits'));
map=fits_read(strcat(savedir,'maps/',dt.name,'_map.fits'));

srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
quad_arr=['A','B','C','D'];
ukmap=zeros(1024);
for iquad=1:4
    quad=quad_arr(iquad);
    sukmaps=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmaps.fits'));
    sukmapg=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapg.fits'));
    sukmap=sukmaps+sukmapg;
    
    if iquad==1
        ukmap(1:512,1:512)=sukmap;
    elseif iquad==2
        ukmap(513:1024,1:512)=sukmap;
    elseif iquad==3
        ukmap(1:512,513:1024)=sukmap;
    else
        ukmap(513:1024,513:1024)=sukmap;
    end
end

%% get diff mask mask to 10 mag
[fr_all] = get_data_frames(inst,dt.name,'flight',flight,'verbose',0);
fr=fr_all(3:end,:,:);
fr1=fr(1:size(fr,1)/2,:,:);
fr2=fr(size(fr,1)/2+1:size(fr,1),:,:);
[map1] = linfit_map(fr1,'verbose',0);
[map2] = linfit_map(fr2,'verbose',0);
diffmap=(map1-map2)./sqrt(2);

strmask=fits_read(strcat(savedir,'masks/',dt.name,'_strmask10.fits')); 
[~,maskin]=get_skymap(diffmap,mask.*strmask,5,3);

strmask=fits_read(strcat(savedir,'masks/',dt.name,'_strmask20.fits')); 
totmask=maskin.*strmask;
clipmax=max(ukmap(find(totmask)));
clipmin=0;

%%
m_arr=10:1:21;

for m_max=m_arr
%%% inst mask * strmask
strmask=fits_read(strcat(savedir,'masks/',dt.name,'_strmask',...
    num2str(m_max),'.fits'));   

totmask=maskin.*strmask;
sp=find(map>clipmax | map<clipmin);
totmask(sp)=0;

mmapcb=(map-mean(map(find(totmask)))).*totmask;
mmapuk=(ukmap-mean(ukmap(find(totmask)))).*totmask;

fits_write(strcat(savedir,'maps/',dt.name,'_cbmap',num2str(m_max),'A'),...
    mmapcb(1:512,1:512));
fits_write(strcat(savedir,'maps/',dt.name,'_cbmap',num2str(m_max),'B'),...
    mmapcb(513:1024,1:512));
fits_write(strcat(savedir,'maps/',dt.name,'_cbmap',num2str(m_max),'C'),...
    mmapcb(1:512,513:1024));
fits_write(strcat(savedir,'maps/',dt.name,'_cbmap',num2str(m_max),'D'),...
    mmapcb(513:1024,513:1024));

fits_write(strcat(savedir,'maps/',dt.name,'_ukmap',num2str(m_max),'A'),...
    mmapuk(1:512,1:512));
fits_write(strcat(savedir,'maps/',dt.name,'_ukmap',num2str(m_max),'B'),...
    mmapuk(513:1024,1:512));
fits_write(strcat(savedir,'maps/',dt.name,'_ukmap',num2str(m_max),'C'),...
    mmapuk(1:512,513:1024));
fits_write(strcat(savedir,'maps/',dt.name,'_ukmap',num2str(m_max),'D'),...
    mmapuk(513:1024,513:1024));

fits_write(strcat(savedir,'masks/',dt.name,'_totmask',num2str(m_max),'A'),...
    totmask(1:512,1:512));
fits_write(strcat(savedir,'masks/',dt.name,'_totmask',num2str(m_max),'B'),...
    totmask(513:1024,1:512));
fits_write(strcat(savedir,'masks/',dt.name,'_totmask',num2str(m_max),'C'),...
    totmask(1:512,513:1024));
fits_write(strcat(savedir,'masks/',dt.name,'_totmask',num2str(m_max),'D'),...
    totmask(513:1024,513:1024));

end


