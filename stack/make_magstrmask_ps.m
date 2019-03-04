flight=40030;
inst=1;
ifield=8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));
m_arr=10:1:20;
%% get inst mask and map (from make_lin_premap.m)
mask=fits_read(strcat(savedir,'masks/',dt.name,'_mask.fits'));
map=fits_read(strcat(savedir,'maps/',dt.name,'_map.fits'));

% correct for the cal factor error
map=map./0.383;
%% masking parameters (mask is saved, don't run it everytime)
%{
for m_max=m_arr

    masks = make_mask_ps(flight,inst,ifield,-1,0,m_max);
    maskg = make_mask_ps(flight,inst,ifield,1,0,m_max);
    strmask = masks.*maskg;

    % Pan-STARRS is incomplete for very bright sources
    tmmask=fits_read(strcat(savedir,'masks2m/',dt.name,'_strmask',...
            num2str(min(m_max,13)),'.fits')); 
    strmask = strmask.*tmmask;

    figure
    imageclip(map.*strmask.*mask);
    title(m_max);

    fits_write(strcat(savedir,'masksps/',dt.name,'_strmask',num2str(m_max)),strmask);
    disp(sprintf('strmask within m=%d save',m_max));
end

%%%get mask of all sources

% the max mag in PanSTARRS is 21.33 (SWIRE_A)
masks = make_mask_ps(flight,inst,ifield,-1,0,30);
maskg = make_mask_ps(flight,inst,ifield,1,0,30);
strmask = masks.*maskg;

tmmask=fits_read(strcat(savedir,'masks2m/',dt.name,'_strmask',...
        num2str(13),'.fits')); 
strmask = strmask.*tmmask;

figure
imageclip(map.*strmask.*mask);

fits_write(strcat(savedir,'masksps/',dt.name,'_strmask_all'),strmask);
disp(sprintf('strmask within m=%d save',m_max));
%}
%% get background level and gradient
strmask = fits_read(strcat(savedir,'masksps/',dt.name,'_strmask_all.fits'));
sigmask = sigclip_mask(map,mask.*strmask,5,5);
grad_map = plane_fit(map,sigmask);
%% get PS sim src map
srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
mapAnames = strcat(srcmapdir,dt.name,'_A_srcmaps_ps.fits');
mapAnameg = strcat(srcmapdir,dt.name,'_A_srcmapg_ps.fits');

psmaps = stick_quad(mapAnames);
psmapg = stick_quad(mapAnameg);
psmap = psmaps + psmapg;
%% subtract mean
totmask=mask.*strmask;
sp=find(map>1.2298e+05 | map<0);
totmask(sp)=0;

cb_bk = median(map(find(totmask)));
ps_bk = median(psmap(find(totmask)));
%%

for m_max=m_arr

strmask=fits_read(strcat(savedir,'masksps/',dt.name,'_strmask',...
    num2str(m_max),'.fits'));   

totmask=mask.*strmask;
mmapcb=(map-grad_map-cb_bk).*totmask;
mmapps=(psmap-ps_bk).*totmask;

figure
setwinsize(gcf,800,300)
subplot(1,2,1)
imageclip(mmapcb);
v = caxis;
subplot(1,2,2)
imageclip(mmapps);
caxis(v);

fits_write(strcat(savedir,'mapsps/',dt.name,'_cbmap',num2str(m_max),'A'),...
    mmapcb(1:512,1:512));
fits_write(strcat(savedir,'mapsps/',dt.name,'_cbmap',num2str(m_max),'B'),...
    mmapcb(513:1024,1:512));
fits_write(strcat(savedir,'mapsps/',dt.name,'_cbmap',num2str(m_max),'C'),...
    mmapcb(1:512,513:1024));
fits_write(strcat(savedir,'mapsps/',dt.name,'_cbmap',num2str(m_max),'D'),...
    mmapcb(513:1024,513:1024));

fits_write(strcat(savedir,'mapsps/',dt.name,'_psmap',num2str(m_max),'A'),...
    mmapps(1:512,1:512));
fits_write(strcat(savedir,'mapsps/',dt.name,'_psmap',num2str(m_max),'B'),...
    mmapps(513:1024,1:512));
fits_write(strcat(savedir,'mapsps/',dt.name,'_psmap',num2str(m_max),'C'),...
    mmapps(1:512,513:1024));
fits_write(strcat(savedir,'mapsps/',dt.name,'_psmap',num2str(m_max),'D'),...
    mmapps(513:1024,513:1024));

fits_write(strcat(savedir,'masksps/',dt.name,'_totmask',num2str(m_max),'A'),...
    totmask(1:512,1:512));
fits_write(strcat(savedir,'masksps/',dt.name,'_totmask',num2str(m_max),'B'),...
    totmask(513:1024,1:512));
fits_write(strcat(savedir,'masksps/',dt.name,'_totmask',num2str(m_max),'C'),...
    totmask(1:512,513:1024));
fits_write(strcat(savedir,'masksps/',dt.name,'_totmask',num2str(m_max),'D'),...
    totmask(513:1024,513:1024));

end


