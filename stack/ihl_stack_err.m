function ihl_stack_err(flight,inst,ifield,varargin)

p = inputParser;

p.addRequired('flight',@isnumeric);
p.addRequired('inst',@isnumeric);
p.addRequired('ifield',@isnumeric);
p.addOptional('addnoise',0,@isnumeric);

p.parse(flight,inst,ifield,varargin{:});

flight = p.Results.flight;
inst = p.Results.inst;
ifield = p.Results.ifield;
addnoise = p.Results.addnoise;


mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
load(sprintf('%s/stackcountdat_%s',loaddir,dt.name),'stackcountdat');

dx = 1200;
verbose = 1;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
mask_inst = stackmapdat(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;
m_min_arr = stackmapdat(ifield).m_min_arr;
m_max_arr = stackmapdat(ifield).m_max_arr;

nbins = 25;

if inst == 1
    stackband = 'I';
else
    stackband = 'H';
end

if addnoise ~= 0
    load(strcat(loaddir,'darkstat'),'darkstat');
    load(strcat(loaddir,'mkkdat'),'mkkdat');
    cp=get_cal_params('flight',flight);
    cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
    frate=cp(inst).framerate;
    slopemean = mean(stackmapdat(ifield).rawmap(find(mkkdat.auto(ifield).mask)));
    slopemap = ones(1024).*slopemean;
    G1 = -3.6;
    pixscale = 7;
    
    phmap=photonnoise_realization(slopemap,G1,dt.nfr,frate);
    rnmap=readnoise_realization(darkstat(ifield).Cl2d_std,pixscale,'norand',1); 
    
    psmap = psmap + phmap.*cal + rnmap.*cal;
    psmap = psmap - mean(psmap(find(mask_inst.*strmask)));
end
    
for im=10:12
m_min = m_min_arr(im);
m_max = m_max_arr(im);

counts = stackcountdat(im).counts;
countg = stackcountdat(im).countg;
dat.counts = counts;
dat.countg = countg;
for type = [1,-1]

if type == 1
    if im ==10
        Nsub = 4;
    else 
        Nsub = 10;
    end
    count_sub = numel(1:Nsub:floor(countg/Nsub)*Nsub);
    dat.countg_sub = count_sub;
    dat.Nsubg = Nsub;
else
    Nsub = 10;
    count_sub = numel(1:Nsub:floor(counts/Nsub)*Nsub);
    dat.counts_sub = count_sub;
    dat.Nsubs = Nsub;
end

[cliplim_arrin,binedgesin,cliplim_arrout,binedgesout]=...
    stackihl_ps_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum,stackband);

profcb_mat = zeros([Nsub, nbins]);
profcberr_mat = zeros([Nsub, nbins]);
profps_mat = zeros([Nsub, nbins]);
profpserr_mat = zeros([Nsub, nbins]);

for isub = 1:Nsub
    [stampercbc,stamperpsc,hitmapc,idx_stack_arr]=...
        stackihl_ps_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,0,verbose,isub,Nsub,stackband);
    
    [stampercb,stamperps,hitmap]=...
        stackihl_ps0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
        mask_inst,strmask,strnum,0,0,verbose,isub,Nsub,stackband,idx_stack_arr);

    stampercb(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stampercbc;
    stamperps(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = stamperpsc;
    hitmap(dx+1-50:dx+1+50,dx+1-50:dx+1+50) = hitmapc;
    
    stackcb = stampercb./hitmap;
    prof = radial_prof(stackcb,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        'weight',hitmap);
    profcb_mat(isub,:) = prof.prof;
    profcberr_mat(isub,:) = prof.err;
    stackps = stamperps./hitmap;
    prof = radial_prof(stackps,ones(2*dx+1),dx+1,dx+1,1,nbins,...
        'weight',hitmap);
    profps_mat(isub,:) = prof.prof;
    profpserr_mat(isub,:) = prof.err;  
end

dat.r_arr = prof.r;
if type == 1
    dat.profgcb_mat = profcb_mat;
    dat.profgps_mat = profps_mat;
    dat.profgcberr_mat = profcberr_mat;
    dat.profgpserr_mat = profpserr_mat;
    dat.profgcb_std = nanstd(profcb_mat);
    dat.profgps_std = nanstd(profps_mat);
else
    dat.profscb_mat = profcb_mat;
    dat.profsps_mat = profps_mat;
    dat.profscberr_mat = profcberr_mat;
    dat.profspserr_mat = profpserr_mat;
    dat.profscb_std = nanstd(profcb_mat);
    dat.profsps_std = nanstd(profps_mat);    
end

end
mcerrdat(im).dat = dat;
end

loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
if addnoise == 0
    save(sprintf('%s/%s_mcerrdat',loaddir,dt.name),'mcerrdat');
else
    save(sprintf('%s/%s_mcerrdat_addnoise',loaddir,dt.name),'mcerrdat');    
end

return