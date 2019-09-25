%%%% run bright source stack for a range of mag bins %%%%%%
function run_psf19_20(inst,ifield)
flight = 40030;
mypaths=get_paths(flight);
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
m_min_arr = [19];
m_max_arr = [20];
%%
dt=get_dark_times(flight,inst,ifield);
cbmap = stackmapdat(ifield).cbmap;
psmap = stackmapdat(ifield).psmap;

im = 1;

m_min = m_min_arr(im);
m_max = m_max_arr(im);
mask_inst = zeros([2,1024,1024]);
mask_inst(1,:,:) = stackmapdat1(ifield).mask_inst_clip;
mask_inst(2,:,:) = stackmapdat2(ifield).mask_inst_clip;
strmask = stackmapdat(ifield).strmask;
strnum = stackmapdat(ifield).strnum;    

psfdat.m_min = m_min;
psfdat.m_max = m_max;

Njack = 16;
srcdat = ps_src_select(flight,inst,ifield,m_min,m_max,mask_inst,...
'sample_type','jack_random','Nsub',Njack);
% only stack 1000 src to speed up
% if srcdat.Ns>1000
%     Nstot = 0;
%     for isub=1:Njack
% %         Nisub = round(srcdat.sub(isub).Ns*1000/srcdat.Ns);
%         Nisub = 2;
%         srcdat.sub(isub).xs_arr = srcdat.sub(isub).xs_arr(1:Nisub);
%         srcdat.sub(isub).ys_arr = srcdat.sub(isub).ys_arr(1:Nisub);
%         srcdat.sub(isub).ms_arr = srcdat.sub(isub).ms_arr(1:Nisub);
%         srcdat.sub(isub).Ns = Nisub;
%         Nstot = Nstot + Nisub;
%     end
%     srcdat.Ns = Nstot;
% end

[clipmaxs, clipmins, r_arr]=...
stackihl_ps0_cliplim(flight,inst,ifield,m_min,m_max,cbmap,psmap,...
mask_inst,strnum,1000,verbose,[],nan,false);

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
psfdat.all.counts = counts;

%%% profile of jackknife samples (leave one out)
for isub=1:Njack
    jackcbs = profcbs - psfdat.sub(isub).profcbs.*psfdat.sub(isub).profhits;
    jackpss = profpss - psfdat.sub(isub).profpss.*psfdat.sub(isub).profhits;
    jackhits = profhits - psfdat.sub(isub).profhits;
    psfdat.jack(isub).profcbs = jackcbs./jackhits;
    psfdat.jack(isub).profpss = jackpss./jackhits;
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

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/%s_psfdat19_20',savedir,dt.name),'psfdat');
return