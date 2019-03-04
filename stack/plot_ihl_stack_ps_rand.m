%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ihl stacking maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
ifield=8;
npix=1200;

dt=get_dark_times(flight,inst,ifield);

savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/'));
%% get the number counts data
ncountfile=strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/ciber_ps/',dt.name,'_stackcounts.txt');
T=readtable(ncountfile);
mbot_arr=T{:,2}';
counts_arr=T{:,5}';
countg_arr=T{:,7}';
countstot_arr=T{:,4}';
countgtot_arr=T{:,6}';
%% plot the random stacking map different source count
loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
num2str(inst),'/'));

figure
setwinsize(gcf,1000,400)

ymin=[];
idx = 0;
count_arr=sort([counts_arr,countg_arr]);
count_arr=count_arr(find(count_arr));
for count=count_arr(13:-3:3)%count_arr(end:-5:1)
    idx = idx + 1;
    stacks=fitsread(strcat(loaddir,'ciber_ps/',dt.name,...
        '_stampers_randFF',num2str(count),'.fits'));
    mstacks=fitsread(strcat(loaddir,'ciber_ps/',dt.name,...
        '_hitmaps_randFF',num2str(count),'.fits'));

    stackscb=stacks./mstacks;
    profile = radial_prof(stackscb,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',0);
    r_arr=profile.r*0.7;
    profs_arr=(profile.prof);
    errs_arr=profile.err;
    ymin=[ymin profs_arr];
    

    subplot(1,2,1)
    errorbar(r_arr,profs_arr,errs_arr,'.-','color',get_color(idx),'DisplayName',...
        strcat(num2str(count),' sources'));hold on
    
    
    subplot(1,2,2)
        loglog(r_arr,profs_arr,'.-');hold on
        errorbar(r_arr,profs_arr,errs_arr,'.-','color',get_color(idx));

drawnow
end


subplot(1,2,1)
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff
plot([4e-1,7e2],[0,0],'k--')


ymin=ymin(find(ymin>0));
ymin=floor(log10(min(ymin)));
ymin=10^ymin;

subplot(1,2,2)
xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')
suptitle('CIBER (w/o FF corretion) random position stacking');
savename=strcat(savedir,dt.name,'_rprof_randomFF');
print(savename,'-dpng');%close

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rand stack same source count different realization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);
mapdir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

% CIBER map (from make_lin_premap.m)
cbmap = fits_read(strcat(mapdir,'maps/',dt.name,'_map.fits'));
% correct for the cal factor error
cbmap = cbmap./0.383;

srcmapdir = strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
mapAnames = strcat(srcmapdir,dt.name,'_A_srcmaps_ps.fits');
mapAnameg = strcat(srcmapdir,dt.name,'_A_srcmapg_ps.fits');

% PS map
psmaps = stick_quad(mapAnames);
psmapg = stick_quad(mapAnameg);
psmap = psmaps + psmapg;

% masks
mask_inst = fits_read(strcat(mapdir,'masks/',dt.name,'_mask.fits'));
strmask = fits_read(strcat(mapdir,'masksps/',dt.name,'_strmask_all.fits'));
totmask = mask_inst.*strmask;

%%% mean & gradient subtraction
sigmask = sigclip_mask(cbmap,totmask,5,5);
grad_cb = plane_fit(cbmap,sigmask);
grad_ps = plane_fit(psmap,sigmask);
cb_bk = mean(cbmap(find(sigmask)));
ps_bk = mean(psmap(find(sigmask)));

cbmap = cbmap - cb_bk - grad_cb;
psmap = psmap - ps_bk - grad_ps;

sig_sp = find((totmask-sigmask)==1);
mask_inst_clip = mask_inst;
mask_inst_clip(sig_sp)=0;
%% stacking with stackihl_ps_randomN.m
figure
setwinsize(gcf,2000,400)

ymin=[];
for i=1:5
[stampercb,stamperps,maskstamper]=...
    stackihl_ps_randomN(flight,inst,ifield,[1000],npix,...
    cbmap,psmap,mask_inst_clip,strmask,savedir);

    stacks=stampercb;
    mstacks=maskstamper;

    stackscb=stacks./mstacks;
    profile = radial_prof(stackscb,ones(2*npix+1),npix+1,npix+1,1,25,...
        'sig',3,'iter_clip',0);
    r_arr=profile.r*0.7;
    profs_arr=(profile.prof);
    errs_arr=profile.err;
    ymin=[ymin profs_arr];
    

    subplot(1,2,1)
    errorbar(r_arr,profs_arr,errs_arr,'.-','color',get_color(i));hold on
    
    
    subplot(1,2,2)
        loglog(r_arr,profs_arr,'.-');hold on
        errorbar(r_arr,profs_arr,errs_arr,'.-','color',get_color(i));

drawnow
end


subplot(1,2,1)
xlabel('arcsec')
ylabel('<I_{stack}>')
xlim([4e-1,7e2])
set(gca, 'XScale', 'log')
plot([4e-1,7e2],[0,0],'k--')


ymin=ymin(find(ymin>0));
ymin=floor(log10(min(ymin)));
ymin=10^ymin;

subplot(1,2,2)
xlim([4e-1,1e3])
xlabel('arcsec')
ylabel('<I_{stack}>')

suptitle('CIBER 2000 random position stacking');
savename=strcat(savedir,dt.name,'_rprof_random1000');
print(savename,'-dpng');%close



