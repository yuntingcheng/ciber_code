function ihl_stack_err_sim(flight,inst,ifield,set,varargin)

p = inputParser;

p.addRequired('flight',@isnumeric);
p.addRequired('inst',@isnumeric);
p.addRequired('ifield',@isnumeric);
p.addRequired('set',@isnumeric);
p.addOptional('addnoise',0,@isnumeric);

p.parse(flight,inst,ifield,set,varargin{:});

flight = p.Results.flight;
inst = p.Results.inst;
ifield = p.Results.ifield;
set = p.Results.set;


mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdatsim',loaddir),'stackmapdatsim');
load(sprintf('%s/stackcountdatsim_%s_%d',loaddir,dt.name,set),'stackcountdat');

dx = 1200;
verbose = 1;

if set==0
    cbmap = stackmapdatsim(ifield).all.cbmap;
    psmap = stackmapdatsim(ifield).all.psmap;
    mask_inst = stackmapdatsim(ifield).all.mask_inst_clip;
    strmask = stackmapdatsim(ifield).all.strmask;
    strnum = stackmapdatsim(ifield).all.strnum;
elseif set==1
    cbmap = stackmapdatsim(ifield).sub.cbmap;
    psmap = stackmapdatsim(ifield).sub.psmap;
    mask_inst = stackmapdatsim(ifield).sub.mask_inst_clip;
    strmask = stackmapdatsim(ifield).sub.strmask;
    strnum = stackmapdatsim(ifield).sub.strnum;
end

nbins = 25;
    
for im=1:3
m_min = im + 15;
m_max = m_min + 1;

counts = stackcountdat(im).counts;
countg = stackcountdat(im).countg;
dat.counts = counts;
dat.countg = countg;
for type = [1,-1]

if type == 1
    if im ==1
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
    stackihl_sim_cliplim(flight,inst,ifield,type,m_min,m_max,...
    cbmap,psmap,mask_inst,strmask,strnum);

profcb_mat = zeros([Nsub, nbins]);
profcberr_mat = zeros([Nsub, nbins]);
profps_mat = zeros([Nsub, nbins]);
profpserr_mat = zeros([Nsub, nbins]);

for isub = 1:Nsub
    [stampercbc,stamperpsc,hitmapc,idx_stack_arr]=...
        stackihl_sim_cent(flight,inst,ifield,type,m_min,m_max,50,cbmap,psmap,...
        mask_inst,strmask,strnum,cliplim_arrin,binedgesin,...
        cliplim_arrout,binedgesout,0,verbose,isub,Nsub);
    
    [stampercb,stamperps,hitmap]=...
        stackihl_sim0(flight,inst,ifield,type,m_min,m_max,dx,cbmap,psmap,...
        mask_inst,strmask,strnum,0,0,verbose,isub,Nsub,idx_stack_arr);

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
save(sprintf('%s/%s_mcerrdatsim%d',loaddir,dt.name,set),'mcerrdat');

return