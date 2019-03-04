%fitparams = get_default_fitparams;
%fitparams.labdata =1;
%fitparams.nskip =0;
%fitparams.nframes = 10000;

indir = '/Users/ytcheng/Documents/MATLAB/CB_darkdata/06-05-2013/';
%%
scan1 = '06-05-2013_13-53-13';
files1 = dir(strcat(indir,'*',scan1,'*'));
nfr = numel(files1)/2;
frames1=zeros(nfr,1024,1024);
for i=1:nfr
    i/nfr  %for print
    infile = strcat(indir,files1(i).name);
    frame = imrotate(fitsread(infile),270);
    frames1(i,:,:)=frame;
end
%%
for cc=1:1024
    for rr=1:1024
        ts = zeros(1,nfr);

        for i=1:nfr
            ts(i)=frames1(i,cc,rr); 
        end
        plot(ts);
    end
end
%%
dark1 = fastlinefit_frin(frames1,0,0,15);
%fastlinefit_frin(frames,quadrant,nskip,nframes)
%%
scan2 = '06-05-2013_14-23-09';
files2 = dir(strcat(indir,'*',scan2,'*'));
nfr = numel(files)/2;
frames2=zeros(nfr,1024,1024);
for i=1:nfr
    i/nfr
    infile = strcat(indir,files2(i).name);
    frame = imrotate(fitsread(infile),270);
    frames2(i,:,:)=frame;
end

ts = zeros(1,nfr);
cc=152;
rr=444;
for i=1:nfr
    ts(i)=frames2(i,cc,rr);
end

%%
dark2 = fastlinefit_frin(frames2,0,0,20);
%%
cp =get_cal_params;
cal = cp(1).apf2eps*110.2;
darknoise = cal*(dark1-dark2)/2;

fitswrite(darknoise,'/home/pkorngut/projects/projects2015/Jul/Jul28_2015_imnoise/im1_20fr_noise_cal.fits');
