%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make the linearized map for stacking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;

mypaths=get_paths(flight);

cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'FFdat'),'FFdat');
load(strcat(loaddir,'maskdat'),'maskdat');
%%
for ifield=5%8:-1:1
dt=get_dark_times(flight,inst,ifield);
[fr] = get_data_frames (inst,dt.name,'flight',flight,'verbose',0);
if ifield==5
    fr=fr(15:end,:,:);
else
    fr=fr(3:end,:,:);
end

%% mask pixel having crazy time stream
dfr=fr(2:end,:,:)-fr(1:end-1,:,:);

fixfr = fr;
shitmask = ones(1024);
shitcut1=1e4;
shitcut2=5e4;
count1=0;count2=0;count3=0;
for cc=1:1024
for rr=1:1024
    
    d1 = dfr(:,cc,rr);
    sp = find(d1 > shitcut2);
    
    if numel(sp)>1
        shitmask(cc,rr)=0;
        count1=count1+1;
    elseif numel(sp)==1
        f1 = fr(:,cc,rr);
        f1(sp+1:end) = f1(sp+1:end) - 2^16;
        fixfr(:,cc,rr)=f1;
        count2=count2+1;
    elseif numel(find(abs(d1)>median(abs(d1))*100))>0
        shitmask(cc,rr)=0;
        count3=count3+1;
    end
end
end
mask=mask_inst.*shitmask;
%% make the unholy flat
FF=FFdat(ifield).FF;
FFmask=ones(1024);FFmask(find(FF==0))=0;
flat = unholy_map(FF,FFmask,200);
%% line fit
[long,off] = linfit_map(fixfr,'verbose',0);
[short] = linfit_map(fixfr(1:4,:,:),'verbose',0);

longcal = cal*long./flat;
shortcal = cal*short./flat;
%% replace saturated pixel with short line fit
qq= squeeze(fixfr(end,:,:));
qqs = qq -off;
sp = find(qqs < -5000);

satmap = zeros(1024);
satmap(sp) = shortcal(sp);

fixmap = longcal;
fixmap(sp) = shortcal(sp);

negmask=ones(1024);
negmask(find(fixmap<0))=0;
%% plot results
%{

cax=([-400,4000]);

cy=842;
cx=137;
dx=30;

%ax = [cx-dx,cx+dx,cy-dx,cy+dx];
ax = [1,1024,1,1024];

subplot(2,2,1)
imageclip(longcal-median(longcal(use)));
caxis(cax)
axis(ax)
subplot(2,2,2)
imageclip(fixmap-median(fixmap(use)));
caxis(cax)
axis(ax)

subplot(2,2,3)
imageclip(satmap);
%caxis([-300,300])
axis(ax)
%}
%% linearized map mean sub
bigmask = maskdat.mask(ifield).bigmask;
use = find(bigmask.*mask);
qq = fixmap-median(fixmap(use));
%%
outdir='/Users/ytcheng/ciber/doc/20170617_Stacking/maps/quadmaps/';

file = strcat(outdir,dt.name,'_inst',num2str(inst),'_map_A.fits');
fits_write(file,qq(1:512,1:512));
file = strcat(outdir,dt.name,'_inst',num2str(inst),'_map_B.fits');
fits_write(file,qq(513:1024,1:512));
file = strcat(outdir,dt.name,'_inst',num2str(inst),'_map_C.fits');
fits_write(file,qq(1:512,513:1024));
file = strcat(outdir,dt.name,'_inst',num2str(inst),'_map_D.fits');
fits_write(file,qq(513:1024,513:1024));

%% mask (inst_mask* bad time stream mask)
outdir='/Users/ytcheng/ciber/doc/20170617_Stacking/maps/quadmasks/';

file = strcat(outdir,dt.name,'_inst',num2str(inst),'_mask_A.fits');
fits_write(file,mask(1:512,1:512));
file = strcat(outdir,dt.name,'_inst',num2str(inst),'_mask_B.fits');
fits_write(file,mask(513:1024,1:512));
file = strcat(outdir,dt.name,'_inst',num2str(inst),'_mask_C.fits');
fits_write(file,mask(1:512,513:1024));
file = strcat(outdir,dt.name,'_inst',num2str(inst),'_mask_D.fits');
fits_write(file,mask(513:1024,513:1024));


disp(sprintf('field%d linearized map and mask written',ifield));
end