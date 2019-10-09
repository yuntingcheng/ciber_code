function run_psf_stack_map(flight,inst,ifield,varargin)
  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addOptional('sample_type','jack_random',@ischar);
  
  p.parse(flight,inst,ifield,varargin{:});

  flight   = p.Results.flight;
  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  sample_type=p.Results.sample_type;
  clear p varargin;

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
Njack = 50;
m_min_arr = [4,4,4,4,12,13,15,16:19];
m_max_arr = [9,10,11,12,13,14,16,17:20];

cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;
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
    
    if m_min<16
        srcdat = tm_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
        'sample_type',sample_type,'Nsub',Njack);

        if srcdat.Ns==0
            continue
        end
        [clipmaxs, clipmins, r_arr]=...
        stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
        mask_inst,strnum,1000,verbose,[],nan,true);
    else
        srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
        'sample_type',sample_type,'Nsub',Njack);
    
        % only stack Ncut src to speed up
%         Ncut = 3000;
%         if srcdat.Ns>Ncut
%             Nstot = 0;
%             for isub=1:Njack
%                 Nisub = round(srcdat.sub(isub).Ns*Ncut/srcdat.Ns);
%                 srcdat.sub(isub).xs_arr = srcdat.sub(isub).xs_arr(1:Nisub);
%                 srcdat.sub(isub).ys_arr = srcdat.sub(isub).ys_arr(1:Nisub);
%                 srcdat.sub(isub).ms_arr = srcdat.sub(isub).ms_arr(1:Nisub);
%                 srcdat.sub(isub).Ns = Nisub;
%                 Nstot = Nstot + Nisub;
%             end
%             srcdat.Ns = Nstot;
%         end
        [clipmaxs, clipmins, r_arr]=...
        stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
        mask_inst,strnum,1000,verbose,[],nan,false);
    end
    
    psfdat.r_arr = r_arr;
    mask_inst = squeeze(mask_inst(inst,:,:));
    
    for isub=1:Njack
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
        psfdat.sub(isub).profcbs = profcbs;
        psfdat.sub(isub).profpss = profpss;        
        psfdat.sub(isub).profhits = profhits;
    end
    
    %%% profile combining all subset
    profcbs = zeros(size(psfdat.r_arr));
    profpss = zeros(size(psfdat.r_arr));
    profhits = zeros(size(psfdat.r_arr));
    counts = 0;
    for isub=1:Njack
        profcbs = profcbs + ...
            psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        profpss = profpss + ...
            psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        profhits = profhits + psfdat.sub(isub).profhits;
        counts = counts + psfdat.sub(isub).counts;
    end
    psfdat.all.profcbs = profcbs./profhits;
    psfdat.all.profpss = profpss./profhits;
    psfdat.all.profhits = profhits;
    psfdat.all.counts = counts;
    
    %%% profile of jackknife samples (leave one out)
    for isub=1:Njack
        jackcbs = profcbs - psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
        jackpss = profpss - psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
        jackhits = profhits - psfdat.sub(isub).profhits;
        psfdat.jack(isub).profcbs = jackcbs./jackhits;
        psfdat.jack(isub).profpss = jackpss./jackhits;
        psfdat.jack(isub).profhits = jackhits;
    end
    
    %%% error bar with jackknife
    errcbs = zeros(size(psfdat.r_arr));
    errpss = zeros(size(psfdat.r_arr));
    for isub=1:Njack
        errcbs = errcbs + ...
            (psfdat.jack(isub).profcbs - psfdat.all.profcbs).^2;        
        errpss = errpss + ...
            (psfdat.jack(isub).profpss - psfdat.all.profpss).^2;
    end
    psfdat.errjack.profcbs = sqrt(errcbs.*((Njack-1)/Njack));
    psfdat.errjack.profpss = sqrt(errpss.*((Njack-1)/Njack));
    psfdatall.mag(im).psfdat = psfdat;
end
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/psfdat_%s',savedir,dt.name),'psfdatall');

%% combine the profile

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/psfdat_%s',savedir,dt.name),'psfdatall');
r_arr = psfdatall.mag(11).psfdat.r_arr;
sp100 = find(r_arr>100);
[rsub_arr,r100] = profile_radial_binning...
    (r_arr,psfdatall.mag(numel(m_min_arr)).psfdat.all.profhits,sp100);

snrs = [];
for im= 1:4
    psfdat = psfdatall.mag(im).psfdat;
    snr = psfdat.all.profcbs./psfdat.errjack.profcbs;
    snrs = [snrs; snr(15:18)];
end
[~,im_out] = max(sum(snrs'));

for im=1:4

im_in = im + 7;
im_mid = 6;
psfdatall.comb(im).r_arr = r_arr;
psfdatall.comb(im).im_in = im_in;
psfdatall.comb(im).im_mid = im_mid;
psfdatall.comb(im).im_out = im_out;

psfdatin = psfdatall.mag(im_in).psfdat;
psfdatmid = psfdatall.mag(im_mid).psfdat;
psfdatout = psfdatall.mag(im_out).psfdat;

%%% combined PSF %%%%
hit = psfdatin.all.profhits;
hit(9:11) = psfdatmid.all.profhits(9:11);
hit(11:end) = psfdatmid.all.profhits(11:end);

psfin = psfdatin.all.profcbs;
psfmid = psfdatmid.all.profcbs;
psfout = psfdatout.all.profcbs;
psfin = psfin./psfin(1);
psfmid = psfmid./psfmid(9).*psfin(9);
psfout = psfout./psfout(11).*psfmid(11);
psfcomb = psfin;
psfcomb(9:11) = psfmid(9:11);
psfcomb(11:end) = psfout(11:end);
[prof15,prof100] = profile_radial_binning(psfcomb,hit,sp100);
psfdatall.comb(im).all.profcb = psfcomb;
psfdatall.comb(im).all.profcbsub = prof15;
psfdatall.comb(im).all.profcb100 = prof100;

psfin = psfdatin.all.profpss;
psfmid = psfdatmid.all.profpss;
psfout = psfdatout.all.profpss;
psfin = psfin./psfin(1);
psfmid = psfmid./psfmid(9).*psfin(9);
psfout = psfout./psfout(11).*psfmid(11);
psfcomb = psfin;
psfcomb(9:11) = psfmid(9:11);
psfcomb(11:end) = psfout(11:end);
[prof15,prof100] = profile_radial_binning(psfcomb,hit,sp100);
psfdatall.comb(im).all.profps = psfcomb;
psfdatall.comb(im).all.profpssub = prof15;
psfdatall.comb(im).all.profps100 = prof100;

psfdatall.comb(im).rsub_arr = rsub_arr;
psfdatall.comb(im).r100 = r100;

%%% combined PSF jackknife%%%%
for isub=1:Njack
    hit = psfdatin.jack(isub).profhits;
    hit(9:11) = psfdatmid.jack(isub).profhits(9:11);
    hit(11:end) = psfdatmid.jack(isub).profhits(11:end);

    psfin = psfdatin.jack(isub).profcbs;
    psfmid = psfdatmid.jack(isub).profcbs;
    psfout = psfdatout.jack(isub).profcbs;
    psfin = psfin./psfin(1);
    psfmid = psfmid./psfmid(9).*psfin(9);
    psfout = psfout./psfout(11).*psfmid(11);
    psfcomb = psfin;
    psfcomb(9:11) = psfmid(9:11);
    psfcomb(11:end) = psfout(11:end);
    [prof15,prof100] = profile_radial_binning(psfcomb,hit,sp100);
    psfdatall.comb(im).jack(isub).profcb = psfcomb;
    psfdatall.comb(im).jack(isub).profcbsub = prof15;
    psfdatall.comb(im).jack(isub).profcb100 = prof100;
    
    psfin = psfdatin.jack(isub).profpss;
    psfmid = psfdatmid.jack(isub).profpss;
    psfout = psfdatout.jack(isub).profpss;
    psfin = psfin./psfin(1);
    psfmid = psfmid./psfmid(9).*psfin(9);
    psfout = psfout./psfout(11).*psfmid(11);
    psfcomb = psfin;
    psfcomb(9:11) = psfmid(9:11);
    psfcomb(11:end) = psfout(11:end);
    [prof15,prof100] = profile_radial_binning(psfcomb,hit,sp100);
    psfdatall.comb(im).jack(isub).profps = psfcomb;
    psfdatall.comb(im).jack(isub).profpssub = prof15;
    psfdatall.comb(im).jack(isub).profps100 = prof100;
    
end

end

save(sprintf('%s/psfdat_%s',savedir,dt.name),'psfdatall');

return