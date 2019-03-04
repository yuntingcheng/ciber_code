fitparams=get_default_fitparams;
fitparams.nskip=6;
fitparams.flight=40030;

alldat(1).name = 'DGL';
alldat(2).name = 'NEP';
alldat(3).name = 'Lockman';
alldat(4).name = 'elat10';
alldat(5).name = 'elat30';
alldat(6).name = 'BootesB';
alldat(7).name = 'BootesA';
alldat(8).name = 'SWIRE';
%% YT -  redo this section from your Dark Current templates
instmask = fitsread_orient('/home/pkorngut/projects/projects2015/Jul/Jul15_2015_imagerstart/inst1masks/inst1_fixedmask.fits');
load('/home/pkorngut/ciberstuff_frombkup/ciberrot_bkup_Jan2014/ciberroot_dev_checkin/lookup/imagers/darks/Averaged_TM1_half1.mat')
t1 = nmdark.map;
load('/home/pkorngut/ciberstuff_frombkup/ciberrot_bkup_Jan2014/ciberroot_dev_checkin/lookup/imagers/darks/Averaged_TM1_half2.mat')
t2 = nmdark.map;
raildark = fliplr((t1 +t2)/2);

%% load the stuff in
indir='/data/CIBER/dr/';

for i=1:8
    inname = strcat(indir,'TM1_',alldat(i).name,...
        '_dr150206.mat');
    load(inname)
    alldat(i).rawmap = data.line.slopedata;
    mask = zeros(1024);
    sp = find(data.mask.mask == 0);
    mask(sp) = 1;
    alldat(i).astrometry = data.astrometry;
    alldat(i).bigmask = mask.*instmask;
    alldat(i).bl_l =data.psf.l;
    alldat(i).bl=data.psf.bl;
 
end
%% subtract dark current
% YT -  redo this section from your Dark Current templates

for i=1:8
    i
    alldat(i).dsub = alldat(i).rawmap - raildark;
    
end

%% Load in all of the ancillary stuff
% YT ignore for now

indir = '/home/pkorngut/ciberstuff_frombkup/ciberrot_bkup_Jan2014/ciberroot_dev_checkin/lookup/imagers/zlgrad/TM1/';
indir2='/home/pkorngut/ciberstuff_frombkup/ciberrot_bkup_Jan2014/ciberroot_dev_checkin/lookup/imagers/dgl/TM1/';
indir3='/home/pkorngut/ciberstuff_frombkup/ciberrot_bkup_Jan2014/ciberroot_dev_checkin/lookup/imagers/masks/starmasks/';
indir4 = '/home/pkorngut/ciberstuff_frombkup/ciberrot_bkup_Jan2014/ciberroot_dev_checkin/lookup/imagers/masks/jackmasks/';

for i=1:8
    i
    inname = strcat(indir,'TM1_40030_zlgrad',alldat(i).name,'.fits');
    alldat(i).zlgrad = fliplr(fitsread_orient(inname));
    inname = strcat(indir2,'TM1_40030_',alldat(i).name,'DGL.fits');
    alldat(i).dglmap = fitsread_orient(inname);
    inname = strcat(indir4,'TM1_40030_',alldat(i).name,'jackmask.fits');
    alldat(i).jackmask = fitsread_orient(inname);

end
%% kill some obvious bright stars
% hand masking of stuff
i=4;
srcmask = mask_source(alldat(i).bigmask,843,135,40);
alldat(i).bigmask =alldat(i).bigmask.*srcmask;
srcmask = mask_source(alldat(i).bigmask,295,210,40);
alldat(i).bigmask =alldat(i).bigmask.*srcmask;

i=5;
srcmask = mask_source(alldat(i).bigmask,438,298,40);
alldat(i).bigmask =alldat(i).bigmask.*srcmask;
srcmask = mask_source(alldat(i).bigmask,771,743,40);
alldat(i).bigmask =alldat(i).bigmask.*srcmask;
srcmask = mask_source(alldat(i).bigmask,162,325,40);
alldat(i).bigmask =alldat(i).bigmask.*srcmask;
%%
for i=1:8
    subplot(3,3,i)
    alldat(i).bigmask = alldat(i).bigmask.*alldat(i).jackmask;
    imageclip(alldat(i).bigmask.*alldat(i).rawmap);
end

%% make the flat
% this can be improved 

goods = [0,0,0,1,1,1,1,1];
stack = zeros(1024);
stackmask = stack;

for i=1:8
    if goods(i)
        use = find(alldat(i).bigmask);
        m = mean(alldat(i).dsub(use));
        mnorm = alldat(i).dsub/m;
        znorm = alldat(i).zlgrad./mean(alldat(i).zlgrad(:)); % including ZL gradients which might not work
        mnorm=mnorm./znorm;
        stack = stack+mnorm.*alldat(i).bigmask;
        stackmask = stackmask + alldat(i).bigmask;
    end
end
flat = stack./stackmask;
good = find(instmask);
flat = flat/nanmean(flat(good));
diff = flat-1.0;
death = find(abs(diff) > .4);
instmask(death)=0;
%%
for i=1:8
    alldat(i).bigmask(death)=0;
end
%%
cp = get_cal_params;
cal = cp(1).apf2eps*(110.2 + 0);
for i=1:8
    alldat(i).calmap = (alldat(i).dsub./flat)*cal;  
    death = find(alldat(i).calmap ~= alldat(i).calmap);
    alldat(i).calmap(death)=0;
    
end
%%
for i=1:8
subplot(3,3,i)
imageclip(alldat(i).calmap.*alldat(i).bigmask);
axis square
title(alldat(i).name);
end

%%
%save('/home/pkorngut/projects/projects2015/Aug/Aug19_2015_newprocess/TM1_processtocal_aug19_2015.mat')
