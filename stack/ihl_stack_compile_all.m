function ihl_stack_compile_all(flight,inst,ifield)

mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);

stackdir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));

pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/plots/PanSTARRS/WISEclass/'));

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

npix = 1200;
nsim = 100;

m_min = ihlprofdat.m_min_arr(10);
m_max = ihlprofdat.m_max_arr(12);
counts = sum(ihlprofdat.counts_arr(10:12));
countg = sum(ihlprofdat.countg_arr(10:12));

%%% stars
profcb = zeros(size(ihlprofdat.r_arr));
profcb_err = zeros(size(ihlprofdat.r_arr));
profps = zeros(size(ihlprofdat.r_arr));
profps_err = zeros(size(ihlprofdat.r_arr));
for iter=1:nsim
    M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
        '_',num2str(counts),'_',num2str(iter),'.csv'));

    profcb = profcb + M(:,2)';
    profps = profps + M(:,4)';
    profcb_err = profcb_err + M(:,2)'.^2;
    profps_err = profps_err + M(:,4)'.^2;
end
profcb_err = profcb_err - profcb;
profps_err = profps_err - profps;   
profcb = profcb./nsim;
profps = profps./nsim;
profcb_err = profcb_err./nsim;
profps_err = profps_err./nsim;

ihlprofdat.all.bk.profscb = profcb;
ihlprofdat.all.bk.profscb_err = profcb_err;
ihlprofdat.all.bk.profsps = profps;
ihlprofdat.all.bk.profsps_err = profps_err;

%%% gals
profcb = zeros(size(ihlprofdat.r_arr));
profcb_err = zeros(size(ihlprofdat.r_arr));
profps = zeros(size(ihlprofdat.r_arr));
profps_err = zeros(size(ihlprofdat.r_arr));
for iter=1:nsim
    M = csvread(strcat(stackdir,'bk_ps/',dt.name,'_randprof',...
        '_',num2str(countg),'_',num2str(iter),'.csv'));

    profcb = profcb + M(:,2)';
    profps = profps + M(:,4)';
    profcb_err = profcb_err + M(:,2)'.^2;
    profps_err = profps_err + M(:,4)'.^2;
end
profcb_err = profcb_err - profcb;
profps_err = profps_err - profps;   
profcb = profcb./nsim;
profps = profps./nsim;
profcb_err = profcb_err./nsim;
profps_err = profps_err./nsim;

ihlprofdat.all.bk.profgcb = profcb;
ihlprofdat.all.bk.profgcb_err = profcb_err;
ihlprofdat.all.bk.profgps = profps;
ihlprofdat.all.bk.profgps_err = profps_err;

%%% get the stacking maps %%%
datascb = zeros(2*npix+1, 2*npix+1);
datagcb = zeros(2*npix+1, 2*npix+1);
datasps = zeros(2*npix+1, 2*npix+1);
datagps = zeros(2*npix+1, 2*npix+1);
mdatas = zeros(2*npix+1, 2*npix+1);
mdatag = zeros(2*npix+1, 2*npix+1);
for m=m_min:m_max-1
    datascb = datascb + fitsread(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_stampers',num2str(m),'_',num2str(m+1),'.fits'));
    datagcb = datagcb + fitsread(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_stamperg',num2str(m),'_',num2str(m+1),'.fits'));                    
    datasps = datasps + fitsread(strcat(stackdir,'/panstarrs_ps/',dt.name,...
        '_stampers',num2str(m),'_',num2str(m+1),'.fits'));
    datagps = datagps + fitsread(strcat(stackdir,'/panstarrs_ps/',dt.name,...
        '_stamperg',num2str(m),'_',num2str(m+1),'.fits'));
    mdat = load(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_hitmaps',num2str(m),'_',num2str(m+1),'.mat'));
    mdatas = mdatas + mdat.hitmap;
    mdat = load(strcat(stackdir,'/ciber_ps/',dt.name,...
        '_hitmapg',num2str(m),'_',num2str(m+1),'.mat'));
    mdatag = mdatag + mdat.hitmap;
end
stackscb=datascb./mdatas;
stackgcb=datagcb./mdatag;
stacksps=datasps./mdatas;
stackgps=datagps./mdatag;
    
 
%%% plot the stacking maps %%%
figure
setwinsize(gcf,800,600)
subplot(2,2,1)
imageclip(stackscb(npix-100:npix+100,npix-100:npix+100));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,2)
imageclip(stackgcb(npix-100:npix+100,npix-100:npix+100));
title(strcat('CIBER',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));
subplot(2,2,3)
imageclip(stacksps(npix-100:npix+100,npix-100:npix+100));
title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(counts),{' '},'stars'));
subplot(2,2,4)
imageclip(stackgps(npix-100:npix+100,npix-100:npix+100));
title(strcat('PanSTARRS',{' '}, 'stack',{' '}, num2str(countg),{' '},'gals'));

suptitle(strcat('PanSTARRS',{' '},num2str(m_min),...
    '<mAB(y band)<',num2str(m_max)));

savename=strcat(pltsavedir,dt.name,'_stackmaps',...
    num2str(m_min),'_',num2str(m_max));
print(savename,'-dpng');%close
    
%%% get the stacking profile %%%
prof = radial_prof(stackscb,ones(2*npix+1),npix+1,npix+1,1,25,...
    'sig',3,'iter_clip',3,'weight',mdatas);
bksub = prof.prof - ihlprofdat.all.bk.profscb;
err = sqrt(prof.err.^2 + ihlprofdat.all.bk.profscb_err.^2);
ihlprofdat.all.data.profscb = bksub;
ihlprofdat.all.data.profscb_err= err;
norm = stackscb(npix+1,npix+1);  
ihlprofdat.all.norm.profscb = bksub./norm;
ihlprofdat.all.norm.profscb_err= err./norm;

prof = radial_prof(stacksps,ones(2*npix+1),npix+1,npix+1,1,25,...
    'sig',3,'iter_clip',3,'weight',mdatas);
bksub = prof.prof - ihlprofdat.all.bk.profsps;
err = sqrt(prof.err.^2 + ihlprofdat.all.bk.profsps_err.^2);
ihlprofdat.all.data.profsps = bksub;
ihlprofdat.all.data.profsps_err= err;
norm = stacksps(npix+1,npix+1);  
ihlprofdat.all.norm.profsps = bksub./norm;
ihlprofdat.all.norm.profsps_err= err./norm;

prof = radial_prof(stackgcb,ones(2*npix+1),npix+1,npix+1,1,25,...
    'sig',3,'iter_clip',3,'weight',mdatag);
bksub = prof.prof - ihlprofdat.all.bk.profgcb;
err = sqrt(prof.err.^2 + ihlprofdat.all.bk.profgcb_err.^2);
ihlprofdat.all.data.profgcb = bksub;
ihlprofdat.all.data.profgcb_err= err;
norm = stackgcb(npix+1,npix+1);  
ihlprofdat.all.norm.profgcb = bksub./norm;
ihlprofdat.all.norm.profgcb_err= err./norm;
    
prof = radial_prof(stackgps,ones(2*npix+1),npix+1,npix+1,1,25,...
    'sig',3,'iter_clip',3,'weight',mdatag);
bksub = prof.prof - ihlprofdat.all.bk.profgps;
err = sqrt(prof.err.^2 + ihlprofdat.all.bk.profgps_err.^2);
ihlprofdat.all.data.profgps = bksub;
ihlprofdat.all.data.profgps_err= err;
norm = stackgps(npix+1,npix+1);  
ihlprofdat.all.norm.profgps = bksub./norm;
ihlprofdat.all.norm.profgps_err= err./norm;
    
%%% get the excess profile %%%
diffcb = ihlprofdat.all.norm.profgcb - ihlprofdat.all.norm.profscb;
diffcb_err = sqrt(ihlprofdat.all.norm.profgcb_err.^2 +...
                  ihlprofdat.all.norm.profscb_err.^2);
diffps = ihlprofdat.all.norm.profgps - ihlprofdat.all.norm.profsps;
diffps_err = sqrt(ihlprofdat.all.norm.profgps_err.^2 +...
                  ihlprofdat.all.norm.profsps_err.^2);
diff = diffcb - diffps;
diff_err = sqrt(diffcb_err.^2 + diffps_err.^2);
ihlprofdat.all.excess.diffcb = diffcb;
ihlprofdat.all.excess.diffcb_err = diffcb_err;
ihlprofdat.all.excess.diffps = diffps;
ihlprofdat.all.excess.diffsp_err = diffps_err;
ihlprofdat.all.excess.diff = diff;
ihlprofdat.all.excess.diff_err = diff_err;

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
save(sprintf('%s/%s_ihlprofdat',loaddir,dt.name),'ihlprofdat');

return