load('/home/pkorngut/projects/projects2015/Jul/Jul30_2015_imager_drcombine/fullpipe_Jul30_dr.mat');
%%
for f=4:8;

fitparams = get_default_fitparams;
fitparams.flight = 40030;
%fr1 = get_data_frames(1,'raildark',fitparams);
fr1 = get_data_frames(1,alldat(f).name,fitparams);
fr1=fr1(2:end,:,:);
s=size(fr1);
%%


dfr = fr1*0;
for i=1:s(1)-1
    dfr(i,:,:)= fr1(i+1,:,:) - fr1(i,:,:);
end

fixfr = fr1;
shitmask = ones(1024);
shitcut1=1e4;
shitcut2=5e4;
rmsmap = zeros(1024);

for cc=1:1024
    cc/1024
    for rr=1:1024
        d1 = dfr(:,cc,rr);
        rmsmap(cc,rr) = std(d1(:));
        sp = find(d1 > 5e4);
        f1 = fr1(:,cc,rr);
        f1(sp+1:end) = f1(sp+1:end) - 2^16;
        fixfr(:,cc,rr)=f1;
        bad = find(abs(d1) > shitcut1 );
        if numel(bad) > 1 ;
            shitmask(cc,rr) =0 ;
        end
    end
end
%%
rmsmask = sigmaclipmask(rmsmap,25,5);
%%

[c,off] = fastlinefit_frin(fixfr,0);
[short] = fastlinefit_frin(fixfr(1:4,:,:),0);
%%
im=instmask;
flatmask = im;
death = find(flat ~= flat);
flat(death) = 0;
flatmask(death) = 0 ;
%flatmask(:,625:645) = 0;

unholyflat = unholy_map(flat,flatmask,200);
flat = unholyflat;
%%
cp = get_cal_params;
cal = cp(1).apf2eps*(110.2 + 0);
calc = cal*c./flat;
calshort = cal*short./flat;

%%

qq= squeeze(fixfr(end,:,:));
qqs = qq -off;
sp = find(qqs < -5000);

satmap = c*0;
satmap(sp) = calshort(sp);

fix = calc;
fix(sp) = calshort(sp);


%%
bm = alldat(f).bigmask;

death = find(fix ~= fix);
%fix(death) =0;
im(death) =0;
bm(death) = 0 ;

death = find(shitmask == 0);
im(death) =0;


use = find(im.*bm);

cax=([-400,4000]);

cy=842;
cx=137;
dx=30;

%ax = [cx-dx,cx+dx,cy-dx,cy+dx];
ax = [1,1024,1,1024];

subplot(2,2,1)
imageclip(calc-median(calc(use)));
caxis(cax)
axis(ax)
subplot(2,2,2)
imageclip(fix-median(fix(use)));
caxis(cax)
axis(ax)

subplot(2,2,3)
imageclip(satmap);
%caxis([-300,300])
axis(ax)
finalmap = fix;

s = finalmap;
death = find(s ~= s);
s(death) = 0;

sm = fillpadsmooth(s,im,1);

subplot(1,1,1)
l=real(log10(sm.*im));
imageclip(l);
caxis([2.5,4.5])
%%
sm2 = sm;
sm2(im == 0) = 0;
%%

outdir='/home/pkorngut/projects/projects2017/Mar/Mar29_2017_allfield_quadmaps/astrinputs/';

file = strcat(outdir,alldat(f).name,'_A.fits');
fits_write(file,sm2(1:512,1:512));m
file = strcat(outdir,alldat(f).name,'_B.fits');
fits_write(file,sm2(513:1024,1:512));
file = strcat(outdir,alldat(f).name,'_C.fits');
fits_write(file,sm2(1:512,513:1024));
file = strcat(outdir,alldat(f).name,'_D.fits');
fits_write(file,sm2(513:1024,513:1024));

%%

qq = fix-median(fix(use));

outdir='/home/pkorngut/projects/projects2017/Mar/Mar29_2017_allfield_quadmaps/quadmaps/';

file = strcat(outdir,alldat(f).name,'_nosm_A.fits');
fits_write(file,qq(1:512,1:512));
file = strcat(outdir,alldat(f).name,'_nosm_B.fits');
fits_write(file,qq(513:1024,1:512));
file = strcat(outdir,alldat(f).name,'_nosm_C.fits');
fits_write(file,qq(1:512,513:1024));
file = strcat(outdir,alldat(f).name,'_nosm_D.fits');
fits_write(file,qq(513:1024,513:1024));

qq = im.*shitmask;

outdir='/home/pkorngut/projects/projects2017/Mar/Mar29_2017_allfield_quadmaps/quadmasks/';

file = strcat(outdir,alldat(f).name,'_mask_A.fits');
fits_write(file,qq(1:512,1:512));
file = strcat(outdir,alldat(f).name,'_mask_B.fits');
fits_write(file,qq(513:1024,1:512));
file = strcat(outdir,alldat(f).name,'_mask_C.fits');
fits_write(file,qq(1:512,513:1024));
file = strcat(outdir,alldat(f).name,'_mask_D.fits');
fits_write(file,qq(513:1024,513:1024));

end