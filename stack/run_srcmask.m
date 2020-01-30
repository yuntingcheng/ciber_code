flight = 40030;
inst = 2;
if inst==1
    Ith=0.5;
else
    Ith=0.5;
end

mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'maskdat'),'maskdat');

for ifield = 4:8
dt=get_dark_times(flight,inst,ifield);

srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
num2str(inst),'/');

map = fits_read(strcat(srcmapdir,dt.name,'_srcmap_ps_all.fits'));

fprintf('make mask inst %d, field %d\n',inst,ifield);

m_max=17;
fprintf('ifield %d m_max = %d\n',ifield,m_max);
[tmmask,tmnum] = make_mask_2m(flight,inst,ifield,0,m_max,Ith,...
    'PSmatch',1,'verbose',false);
[psmask,psnum] = make_mask_ps(flight,inst,ifield,0,0,m_max,Ith,...
    'verbose',false);
strmask = psmask.*tmmask;
strnum = psnum + tmnum;
maskdat.mask(ifield).m_max(m_max).strmask_stack = strmask;
maskdat.mask(ifield).m_max(m_max).strnum_stack = strnum;    

for m_max=[18,19,20]
    strmask1 = maskdat.mask(ifield).m_max(m_max-1).strmask_stack;
    strnum1 = maskdat.mask(ifield).m_max(m_max-1).strnum_stack;
    fprintf('ifield %d m_max = %d\n',ifield,m_max);
    [tmmask,tmnum] = make_mask_2m(flight,inst,ifield,m_max-1,m_max,Ith,...
        'PSmatch',1,'verbose',false);
    [psmask,psnum] = make_mask_ps(flight,inst,ifield,0,m_max-1,m_max,Ith,...
        'verbose',false);
    strmask = psmask.*tmmask.*strmask1;
    strnum = psnum + tmnum + strnum1;
    maskdat.mask(ifield).m_max(m_max).strmask_stack = strmask;
    maskdat.mask(ifield).m_max(m_max).strnum_stack = strnum;    
end


fprintf('ifield %d all\n',ifield);
strmask1 = maskdat.mask(ifield).m_max(20).strmask_stack;
strnum1 = maskdat.mask(ifield).m_max(20).strnum_stack;
[tmmask,tmnum] = make_mask_2m(flight,inst,ifield,20,23,Ith,...
    'PSmatch',1,'verbose',false);
[psmask,psnum] = make_mask_ps(flight,inst,ifield,0,20,23,Ith,'verbose',false);

strmask = psmask.*tmmask.*strmask1;
strnum = psnum + tmnum + strnum1;

maskdat.mask(ifield).strmask_stack = strmask;
maskdat.mask(ifield).strnum_stack = strnum;


save(strcat(loaddir,'maskdat'),'maskdat');
end

%% HSC

flight = 40030;
inst = 2;
ifield = 8;
if inst==1
    Ith=0.5;
else
    Ith=0.5;
end
mypaths=get_paths(flight);

for hsc_idx= 0:11
[field,mask_inst] = HSC_fields_info(hsc_idx);
fprintf('make mask inst %d, %s\n',inst,field);
srcmapdir = strcat(mypaths.ciberdir,'doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
srcmap1 = fits_read(strcat(srcmapdir,'srcmap_sim1_all_hsc_',field,'.fits'));
srcmap2 = fits_read(strcat(srcmapdir,'srcmap_sim2_all_hsc_',field,'.fits'));

m_max=17;
fprintf('%s m_max = %d\n',field,m_max);
[strmask,strnum] = make_mask_hsc(flight,inst,ifield,field,-5,m_max,Ith);
maskdat.mask(hsc_idx+1).name = field;
maskdat.mask(hsc_idx+1).m_max(m_max).strmask_stack = strmask;
maskdat.mask(hsc_idx+1).m_max(m_max).strnum_stack = strnum;    

for m_max=[18,19,20]
    strmask1 = maskdat.mask(hsc_idx+1).m_max(m_max-1).strmask_stack;
    strnum1 = maskdat.mask(hsc_idx+1).m_max(m_max-1).strnum_stack;
    fprintf('%s m_max = %d\n',field,m_max);
    [strmask,strnum] = make_mask_hsc(flight,inst,ifield,field,m_max-1,m_max,Ith);
    strmask = strmask.*strmask1;
    strnum = strnum + strnum1;
    maskdat.mask(hsc_idx+1).m_max(m_max).strmask_stack = strmask;
    maskdat.mask(hsc_idx+1).m_max(m_max).strnum_stack = strnum;    
end

fprintf('%s m_max = %d\n',field,20);
strmask1 = maskdat.mask(hsc_idx+1).m_max(20).strmask_stack;
strnum1 = maskdat.mask(hsc_idx+1).m_max(20).strnum_stack;
[strmask,strnum] = make_mask_hsc(flight,inst,ifield,field,20,22,Ith);

strmask = strmask.*strmask1;
strnum = strnum + strnum1;

maskdat.mask(hsc_idx+1).strmask_stack = strmask;
maskdat.mask(hsc_idx+1).strnum_stack = strnum;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(strcat(savedir,'maskdathsc'),'maskdat');
end

% flight = 40030;
% inst = 2;
% ifield = 8;
% mypaths=get_paths(flight);
% savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
% for hsc_idx= 0:11
%     load(strcat(savedir,'maskdathsc',num2str(hsc_idx)),'maskdat');
% 
%     maskdatall.mask(hsc_idx+1) = maskdat.mask(hsc_idx+1);
% end
% maskdat = maskdatall;
% save(strcat(savedir,'maskdathsc'),'maskdat');
