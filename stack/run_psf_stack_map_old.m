function run_psf_stack_map_old(flight,inst,ifield)
%%%%%%%%%%%%%%%%%%%%%%
% stack PSF
% use 2MASS j band, 4<m<12 for wing, and 15<m<16 for core.
% Two profile matched at 11th r_arr bin (r=16.0826'')
%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  
  p.parse(flight,inst,ifield);

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  
  clear p varargin;
%%
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
load(sprintf('%s/TM%d/stackmapdat',mypaths.alldat,1),'stackmapdat');
stackmapdat1 = stackmapdat;
load(sprintf('%s/TM%d/stackmapdat',mypaths.alldat,2),'stackmapdat');
stackmapdat2 = stackmapdat;

if inst==1
    stackmapdat = stackmapdat1;
else
    stackmapdat = stackmapdat2;
end
    
dx = 1200;
verbose = false;
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
m_min_arr = [4,15];
m_max_arr = [12,16];

%%
for im= 1:numel(m_min_arr)

    m_min = m_min_arr(im);
    m_max = m_max_arr(im);
    mask_inst = zeros([2,1024,1024]);
    mask_inst(1,:,:) = stackmapdat1(ifield).mask_inst_clip;
    mask_inst(2,:,:) = stackmapdat2(ifield).mask_inst_clip;
    strmask = stackmapdat(ifield).strmask;
    strnum = stackmapdat(ifield).strnum;    
    
    psfdat.m_min = m_min;
    psfdat.m_max = m_max;
    
    srcdat = tm_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
    'sample_type','jack_random','Nsub',10);
    
    if srcdat.Ns==0
        continue
    end

    [clipmaxs, clipmins, r_arr]=...
    stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
    mask_inst,strnum,1000,verbose,[],nan,true);
    psfdat.r_arr = r_arr;
    mask_inst = squeeze(mask_inst(inst,:,:));
    
    cbmean = mean(cbmap(find(mask_inst.*strmask)));
    psmean = mean(psmap(find(mask_inst.*strmask)));

    for isub=1:10
        [~,~,~,profcbs,profpss,profhits] = ...
            stackihl_ps0_hist_map(flight,inst,ifield,dx,cbmap,psmap,mask_inst,...
            strmask,strnum,1,verbose,nan,clipmaxs,clipmins,...
            srcdat.sub(isub).xs_arr,srcdat.sub(isub).ys_arr,...
            srcdat.sub(isub).ms_arr,true);

        fprintf('stack %s, %d<m<%d, isub %d, %d srcs\n',...
            dt.name,m_min,m_max,isub,srcdat.sub(isub).Ns);

        psfdat.sub(isub).counts = srcdat.sub(isub).Ns;
        profcbs(profhits==0) = 0;
        profpss(profhits==0) = 0;
        psfdat.sub(isub).profcbs = profcbs - cbmean;
        psfdat.sub(isub).profpss = profpss - psmean;        
        psfdat.sub(isub).profhits = profhits;
    end
    
    %%% profile combining all subset
    profcbs = zeros(size(psfdat.r_arr));
    profpss = zeros(size(psfdat.r_arr));
    profhits = zeros(size(psfdat.r_arr));
    counts = 0;
    for isub=1:10
        profcbs = profcbs + ...
            psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        profpss = profpss + ...
            psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        profhits = profhits + psfdat.sub(isub).profhits;
        counts = counts + psfdat.sub(isub).counts;
    end
    psfdat.all.profcbs = profcbs./profhits;
    psfdat.all.profpss = profpss./profhits;
    psfdat.all.counts = counts;
    
    %%% profile of jackknife samples (leave one out)
    for isub=1:10
        jackcbs = profcbs - psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        jackpss = profpss - psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        jackhits = profhits - psfdat.sub(isub).profhits;
        psfdat.jack(isub).profcbs = jackcbs./jackhits;
        psfdat.jack(isub).profpss = jackpss./jackhits;
    end
    
    %%% error bar with jackknife
    errcbs = zeros(size(psfdat.r_arr));
    errpss = zeros(size(psfdat.r_arr));
    for isub=1:10
        errcbs = errcbs + ...
            (psfdat.jack(isub).profcbs - psfdat.all.profcbs).^2;        
        errpss = errpss + ...
            (psfdat.jack(isub).profpss - psfdat.all.profpss).^2;
    end
    psfdat.errjack.profcbs = sqrt(errcbs.*(9/10));
    psfdat.errjack.profpss = sqrt(errpss.*(9/10));
    
    %%% 
    if im==1
        psfcbout = psfdat.all.profcbs./psfdat.all.profcbs(11);
        psfcbout_err = psfdat.errjack.profcbs./psfdat.all.profcbs(11);
        psfpsout = psfdat.all.profpss./psfdat.all.profpss(11);
        psfpsout_err = psfdat.errjack.profpss./psfdat.all.profpss(11);
    else
        psfcbin = psfdat.all.profcbs./psfdat.all.profcbs(11);
        psfcbin_err = psfdat.errjack.profcbs./psfdat.all.profcbs(11);        
        psfpsin = psfdat.all.profpss./psfdat.all.profpss(11);
        psfpsin_err = psfdat.errjack.profpss./psfdat.all.profpss(11);        
    end  
    
end
psfcb = psfcbout;
psfcb_err = psfcbout_err;
psfcb(1:11) = psfcbin(1:11);
psfcb_err(1:11) = psfcbin_err(1:11);
psfcb_err = psfcb_err./psfcb(1);
psfcb = psfcb./psfcb(1);

psfps = psfpsout;
psfps_err = psfpsout_err;
psfps(1:11) = psfpsin(1:11);
psfps_err(1:11) = psfpsin_err(1:11);
psfps_err = psfps_err./psfps(1);
psfps = psfps./psfps(1);
%%
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

for masklim = [0, 1]
    if masklim==0
        load(sprintf('%s/stackdat_%s',savedir,dt.name),'stackdatall');
    else
        load(sprintf('%s/stackdat_%s_masklim',savedir,dt.name),'stackdatall');
    end
    
for im=1:4
    stackdat = stackdatall(im).stackdat;
    stackdat.psfcb = psfcb;
    stackdat.psfcb_err = psfcb_err;
    
    norm = stackdat.norm.profcbg(1);
    stackdat.norm.profcbpsf = psfcb.*norm;
    stackdat.norm.profcbpsf_err = psfcb_err.*norm;
    stackdat.norm.normcbpsf = norm;
    stackdat.norm.profcbpsf100_upperbound = ...
        stackdat.norm.profcbpsf(18)+stackdat.norm.profcbpsf_err(18);
    diffcb = stackdat.norm.profcbg - stackdat.norm.profcbpsf;
    diffcb100 = stackdat.norm.profcbg100;
    diffcb_err = stackdat.norm.profcbg_err;
    diffcb_err100 = stackdat.norm.profcbg_err100;
    stackdat.excess.diffcbpsf = diffcb;
    stackdat.excess.diffcbpsf100 = diffcb100;
    stackdat.excess.diffcbpsf_err = diffcb_err;
    stackdat.excess.diffcbpsf_err100 = diffcb_err100;
    
    norm = stackdat.norm.profpsg(1);
    stackdat.norm.profpspsf = psfps.*norm;
    stackdat.norm.profpspsf_err = psfps_err.*norm;
    stackdat.norm.normpspsf = norm;
    stackdat.norm.profpspsf100_upperbound = ...
        stackdat.norm.profpspsf(18)+stackdat.norm.profpspsf_err(18);
    diffps = stackdat.norm.profpsg - stackdat.norm.profpspsf;
    diffps100 = stackdat.norm.profpsg100;
    diffps_err = stackdat.norm.profpsg_err;
    diffps_err100 = stackdat.norm.profpsg_err100;
    stackdat.excess.diffpspsf = diffps;
    stackdat.excess.diffpspsf100 = diffps100;
    stackdat.excess.diffpspsf_err = diffps_err;
    stackdat.excess.diffpspsf_err100 = diffps_err100;
    
    %%% binning inner and outer bins for cov %%%
    Njack = numel(stackdat.jack);
    w = 0;
    for i=1:Njack
        w = w + stackdat.sub(i).profhitg;
    end
    stackdat.radweight = w;
    rsub_arr = zeros([1,15]);
    rsub_arr(2:end-1) = r_arr(7:19);
    rsub_arr(1) = sum(r_arr(1:6).*w(1:6))./sum(w(1:6));
    rsub_arr(end) = sum(r_arr(20:25).*w(20:25))./sum(w(20:25));
    stackdat.rsub_arr = rsub_arr;
    
    excess_jack = zeros([Njack,numel(r_arr)]);
    excessin = zeros([1,Njack]);
    excessout = zeros([1,Njack]);
    for i=1:Njack
        excess_jack(i,:) = (stackdat.jack(i).profcbg -...
            stackdat.bk.profcbg) - stackdat.norm.profcbpsf;
        excessin(i) = sum(excess_jack(i,1:6).*w(1:6))./sum(w(1:6));
        excessout(i) = sum(excess_jack(i,20:25).*w(20:25))./sum(w(20:25));
    end
    excess_jacksub = zeros([Njack,15]);
    excess_jacksub(:,2:end-1) = excess_jack(:,7:19);
    excess_jacksub(:,1) = excessin;
    excess_jacksub(:,end) = excessout;
    
    %%% get cov %%%
    cov_mat = zeros(numel(r_arr));
    cov_rho = zeros(numel(r_arr));
    for i=1:numel(r_arr)
        for j=i:numel(r_arr)
            dati = excess_jack(:,i);
            datj = excess_jack(:,j);
            cov_mat(i,j) = mean(dati.*datj) - mean(dati)*mean(datj);
            cov_rho(i,j) = cov_mat(i,j)./sqrt(var(dati)*var(datj));
            cov_mat(j,i) = cov_mat(i,j);
            cov_rho(j,i) = cov_rho(i,j);
        end
    end
    cov_mat = cov_mat.*(Njack - 1);
    cov_inv = inv(cov_mat);
    stackdat.cov.cov_mat = cov_mat;
    stackdat.cov.cov_rho = cov_rho;
    stackdat.cov.cov_inv = cov_inv;
    
    cov_mat = zeros(numel(rsub_arr));
    cov_rho = zeros(numel(rsub_arr));
    for i=1:numel(rsub_arr)
        for j=i:numel(rsub_arr)
            dati = excess_jacksub(:,i);
            datj = excess_jacksub(:,j);
            cov_mat(i,j) = mean(dati.*datj) - mean(dati)*mean(datj);
            cov_rho(i,j) = cov_mat(i,j)./sqrt(var(dati)*var(datj));
            cov_mat(j,i) = cov_mat(i,j);
            cov_rho(j,i) = cov_rho(i,j);
        end
    end
    cov_mat = cov_mat.*(Njack - 1);
    cov_inv = inv(cov_mat);
    stackdat.cov.cov_matsub = cov_mat;
    stackdat.cov.cov_rhosub = cov_rho;
    stackdat.cov.cov_invsub = cov_inv;
     
    stackdatall(im).stackdat = stackdat;
end

if masklim==0
    save(sprintf('%s/stackdat_%s',savedir,dt.name),'stackdatall');
else
    save(sprintf('%s/stackdat_%s_masklim',savedir,dt.name),'stackdatall');
end
    
end
return