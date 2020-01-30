function run_excess_stack_map(flight,inst,ifield,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addOptional('masklim',false,@islogical);
  
  p.parse(flight,inst,ifield,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  masklim = p.Results.masklim;
  clear p varargin;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
if masklim
    load(sprintf('%s/stackdat_%s_masklim',savedir,dt.name),'stackdatall');        
else
    load(sprintf('%s/stackdat_%s',savedir,dt.name),'stackdatall');
end
load(sprintf('%s/psfdat_%s',savedir,dt.name),'psfdatall');

for im=1:4
psfdat = psfdatall.comb(im);
stackdat = stackdatall(im).stackdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
r_arr = stackdat.r_arr;
rsub_arr = stackdat.rsub_arr;
Njack = numel(stackdat.sub);

excessdat.m_min = m_min;
excessdat.m_max = m_max;
excessdat.r_arr = r_arr;
excessdat.rsub_arr = rsub_arr;
excessdat.r100 = stackdat.r100;
excessdat.dat = stackdat.all;
excessdat.bg = stackdat.bg;
excessdat.psf = psfdat.all;

normcb = excessdat.dat.profcbg(1) - excessdat.bg.profcbg(1);
normps = excessdat.dat.profpsg(1) - excessdat.bg.profpsg(1);
excessdat.normcb = normcb;
excessdat.normps = normps;

%%%% excess from individual terms %%%%%
excessdat.excess.profcbg = (excessdat.dat.profcbg - ...
    excessdat.bg.profcbg) - normcb.*excessdat.psf.profcb;
excessdat.excess.profpsg = (excessdat.dat.profpsg - ...
    excessdat.bg.profpsg) - normps.*excessdat.psf.profps;
excessdat.excess.profcbgsub = (excessdat.dat.profcbgsub - ...
    excessdat.bg.profcbgsub) - normcb.*excessdat.psf.profcbsub;
excessdat.excess.profpsgsub = (excessdat.dat.profpsgsub - ...
    excessdat.bg.profpsgsub) - normps.*excessdat.psf.profpssub;
excessdat.excess.profcbg100 = (excessdat.dat.profcbg100 - ...
    excessdat.bg.profcbg100) - normcb.*excessdat.psf.profcb100;
excessdat.excess.profpslg100 = (excessdat.dat.profpsg100 - ...
    excessdat.bg.profpsg100) - normps.*excessdat.psf.profps100;

%%%%%% jackknife of individual terms %%%%%%%
dat_profcb = zeros([Njack,numel(r_arr)]);
dat_profps = zeros([Njack,numel(r_arr)]);
dat_profcbsub = zeros([Njack,numel(rsub_arr)]);
dat_profpssub = zeros([Njack,numel(rsub_arr)]);
dat_profcb100 = zeros([Njack,1]);
dat_profps100 = zeros([Njack,1]); 
for isub=1:Njack
    dat_profcb(isub,:) = stackdat.jack(isub).profcbg;
    dat_profps(isub,:) = stackdat.jack(isub).profpsg;
    dat_profcbsub(isub,:) = stackdat.jack(isub).profcbgsub;
    dat_profpssub(isub,:) = stackdat.jack(isub).profpsgsub;
    dat_profcb100(isub) = stackdat.jack(isub).profcbg100;
    dat_profps100(isub) = stackdat.jack(isub).profpsg100;
end
excessdat.datcov.covcb = get_cov_matrix(dat_profcb).*(Njack-1);
excessdat.datcov.covcbsub = get_cov_matrix(dat_profcbsub).*(Njack-1);
excessdat.datcov.covcb100 = get_cov_matrix(dat_profcb100).*(Njack-1);
excessdat.datcov.covps = get_cov_matrix(dat_profps).*(Njack-1);
excessdat.datcov.covpssub = get_cov_matrix(dat_profpssub).*(Njack-1);
excessdat.datcov.covps100 = get_cov_matrix(dat_profps100).*(Njack-1);

dat_profcb = zeros([Njack,numel(r_arr)]);
dat_profps = zeros([Njack,numel(r_arr)]);
dat_profcbsub = zeros([Njack,numel(rsub_arr)]);
dat_profpssub = zeros([Njack,numel(rsub_arr)]);
dat_profcb100 = zeros([Njack,1]);
dat_profps100 = zeros([Njack,1]); 
for isub=1:Njack
    dat_profcb(isub,:) = psfdat.jack(isub).profcb.*normcb;
    dat_profps(isub,:) = psfdat.jack(isub).profps.*normps;
    dat_profcbsub(isub,:) = psfdat.jack(isub).profcbsub.*normcb;
    dat_profpssub(isub,:) = psfdat.jack(isub).profpssub.*normps;
    dat_profcb100(isub) = psfdat.jack(isub).profcb100.*normcb;
    dat_profps100(isub) = psfdat.jack(isub).profps100.*normps;
end
covcb = get_cov_matrix(dat_profcb).*(Njack-1);
covcb(1,:) = 0;
covcb(:,1) = 0;
excessdat.psfcov.covcb = covcb;
excessdat.psfcov.covcbsub = get_cov_matrix(dat_profcbsub).*(Njack-1);
excessdat.psfcov.covcb100 = get_cov_matrix(dat_profcb100).*(Njack-1);
covps = get_cov_matrix(dat_profps).*(Njack-1);
covps(1,:) = 0;
covps(:,1) = 0;
excessdat.psfcov.covps = covps;
excessdat.psfcov.covpssub = get_cov_matrix(dat_profpssub).*(Njack-1);
excessdat.psfcov.covps100 = get_cov_matrix(dat_profps100).*(Njack-1);

excessdat.bgcov.covcb = stackdat.bgcov.covcb;
excessdat.bgcov.covcbsub = stackdat.bgcov.covcbsub;
excessdat.bgcov.covcb100 = stackdat.bgcov.covcb100;
excessdat.bgcov.covps = stackdat.bgcov.covps;
excessdat.bgcov.covpssub = stackdat.bgcov.covpssub;
excessdat.bgcov.covps100 = stackdat.bgcov.covps100;

excessdat.excov.covcb = excessdat.datcov.covcb + ...
    excessdat.bgcov.covcb + excessdat.psfcov.covcb;
excessdat.excov.covcbsub = excessdat.datcov.covcbsub + ...
    excessdat.bgcov.covcbsub + excessdat.psfcov.covcbsub;
excessdat.excov.covcb100 = excessdat.datcov.covcb100 + ...
    excessdat.bgcov.covcb100 + excessdat.psfcov.covcb100;
excessdat.excov.covps = excessdat.datcov.covps + ...
    excessdat.bgcov.covps + excessdat.psfcov.covps;
excessdat.excov.covpssub = excessdat.datcov.covpssub + ...
    excessdat.bgcov.covpssub + excessdat.psfcov.covpssub;
excessdat.excov.covps100 = excessdat.datcov.covps100 + ...
    excessdat.bgcov.covps100 + excessdat.psfcov.covps100;

%%%% excess from jackknife excess w/ same PSF norm%%%%%
dat_profcb = zeros([Njack,numel(r_arr)]);
dat_profps = zeros([Njack,numel(r_arr)]);
dat_profcbsub = zeros([Njack,numel(rsub_arr)]);
dat_profpssub = zeros([Njack,numel(rsub_arr)]);
dat_profcb100 = zeros([Njack,1]);
dat_profps100 = zeros([Njack,1]);

bgcb = excessdat.bg.profcbg;
bgcbsub = excessdat.bg.profcbgsub;
bgcb100 = excessdat.bg.profcbg100;
bgps = excessdat.bg.profpsg;
bgpssub = excessdat.bg.profpsgsub;
bgps100 = excessdat.bg.profpsg100;
for isub=1:Njack
    datcb = stackdat.jack(isub).profcbg;
    psfcb = psfdat.jack(isub).profcb.*normcb;
    datcbsub = stackdat.jack(isub).profcbgsub;
    psfcbsub = psfdat.jack(isub).profcbsub.*normcb;
    datcb100 = stackdat.jack(isub).profcbg100;
    psfcb100 = psfdat.jack(isub).profcb100.*normcb;
    
    datps = stackdat.jack(isub).profpsg;
    psfps = psfdat.jack(isub).profps.*normps;
    datpssub = stackdat.jack(isub).profpsgsub;
    psfpssub = psfdat.jack(isub).profpssub.*normps;
    datps100 = stackdat.jack(isub).profpsg100;
    psfps100 = psfdat.jack(isub).profps100.*normps;
    
    dat_profcb(isub,:) = datcb - bgcb - psfcb;
    dat_profps(isub,:) = datps - bgps - psfps;
    dat_profcbsub(isub,:) = datcbsub - bgcbsub - psfcbsub;
    dat_profpssub(isub,:) = datpssub - bgpssub - psfpssub;
    dat_profcb100(isub) = datcb100 - bgcb100 - psfcb100;
    dat_profps100(isub) = datps100 - bgps100 - psfps100;
end
excessdat.exjcov.covcb = excessdat.bgcov.covcb + ...
    get_cov_matrix(dat_profcb).*(Njack-1);
excessdat.exjcov.covcbsub = excessdat.bgcov.covcbsub + ...
    get_cov_matrix(dat_profcbsub).*(Njack-1);
excessdat.exjcov.covcb100 = excessdat.bgcov.covcb100 + ...
    get_cov_matrix(dat_profcb100).*(Njack-1);
excessdat.exjcov.covps = excessdat.bgcov.covps + ...
    get_cov_matrix(dat_profps).*(Njack-1);
excessdat.exjcov.covpssub = excessdat.bgcov.covpssub + ...
    get_cov_matrix(dat_profpssub).*(Njack-1);
excessdat.exjcov.covps100 = excessdat.bgcov.covps100 + ...
    get_cov_matrix(dat_profps100).*(Njack-1);

excessdat.excessj.profcbg = mean(dat_profcb);
excessdat.excessj.profpsg = mean(dat_profps);
excessdat.excessj.profcbgsub = mean(dat_profcbsub);
excessdat.excessj.profpsgsub = mean(dat_profpssub);
excessdat.excessj.profcbg100 = mean(dat_profcb100);
excessdat.excessj.profpsg100 = mean(dat_profps100);

%%%% excess from jackknife excess w/ same PSF norm%%%%%
dat_profcb = zeros([Njack,numel(r_arr)]);
dat_profps = zeros([Njack,numel(r_arr)]);
dat_profcbsub = zeros([Njack,numel(rsub_arr)]);
dat_profpssub = zeros([Njack,numel(rsub_arr)]);
dat_profcb100 = zeros([Njack,1]);
dat_profps100 = zeros([Njack,1]); 

bgcb = excessdat.bg.profcbg;
bgcbsub = excessdat.bg.profcbgsub;
bgcb100 = excessdat.bg.profcbg100;
bgps = excessdat.bg.profpsg;
bgpssub = excessdat.bg.profpsgsub;
bgps100 = excessdat.bg.profpsg100;
for isub=1:Njack
    datcb = stackdat.jack(isub).profcbg;
    normcb = datcb(1) - bgcb(1);
    psfcb = psfdat.jack(isub).profcb.*normcb;
    datcbsub = stackdat.jack(isub).profcbgsub;
    psfcbsub = psfdat.jack(isub).profcbsub.*normcb;
    datcb100 = stackdat.jack(isub).profcbg100;
    psfcb100 = psfdat.jack(isub).profcb100.*normcb;
    
    datps = stackdat.jack(isub).profpsg;
    normps = datps(1) - bgps(1);
    psfps = psfdat.jack(isub).profps.*normps;
    datpssub = stackdat.jack(isub).profpsgsub;
    psfpssub = psfdat.jack(isub).profpssub.*normps;
    datps100 = stackdat.jack(isub).profpsg100;
    psfps100 = psfdat.jack(isub).profps100.*normps;
    
    dat_profcb(isub,:) = datcb - bgcb - psfcb;
    dat_profps(isub,:) = datps - bgps - psfps;
    dat_profcbsub(isub,:) = datcbsub - bgcbsub - psfcbsub;
    dat_profpssub(isub,:) = datpssub - bgpssub - psfpssub;
    dat_profcb100(isub) = datcb100 - bgcb100 - psfcb100;
    dat_profps100(isub) = datps100 - bgps100 - psfps100;
end
excessdat.exjicov.covcb = excessdat.bgcov.covcb + ...
    get_cov_matrix(dat_profcb).*(Njack-1);
excessdat.exjicov.covcbsub = excessdat.bgcov.covcbsub + ...
    get_cov_matrix(dat_profcbsub).*(Njack-1);
excessdat.exjicov.covcb100 = excessdat.bgcov.covcb100 + ...
    get_cov_matrix(dat_profcb100).*(Njack-1);
excessdat.exjicov.covps = excessdat.bgcov.covps + ...
    get_cov_matrix(dat_profps).*(Njack-1);
excessdat.exjicov.covpssub = excessdat.bgcov.covpssub + ...
    get_cov_matrix(dat_profpssub).*(Njack-1);
excessdat.exjicov.covps100 = excessdat.bgcov.covps100 + ...
    get_cov_matrix(dat_profps100).*(Njack-1);

excessdat.excessji.profcbg = mean(dat_profcb);
excessdat.excessji.profpsg = mean(dat_profps);
excessdat.excessji.profcbgsub = mean(dat_profcbsub);
excessdat.excessji.profpsgsub = mean(dat_profpssub);
excessdat.excessji.profcbg100 = mean(dat_profcb100);
excessdat.excessji.profpsg100 = mean(dat_profps100);

excessdatall(im).excessdat = excessdat;
end

if masklim
    save(sprintf('%s/excessdat_%s_masklim',...
        savedir,dt.name),'excessdatall');        
else
    save(sprintf('%s/excessdat_%s',...
        savedir,dt.name),'excessdatall');
end
clear excessdat

return