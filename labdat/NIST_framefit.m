%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find reset of NIST imager data and save.
% Save the reset info, frames to used, line fit map in struc.
% modified from find_reset.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% I band

time1={'10-40-20';'10-41-15'};
time2={'10-42-59';'10-43-23';'10-44-52';'10-45-13';'10-45-40'};
time3={'10-50-45';'10-51-10';'10-51-38';'10-52-12';'10-52-37'};
time4={'10-54-11';'10-55-10'};
time_arr={time1,time2,time3,time4};

loaddir='/Volumes/HD1TB/CIBER/data/Cal_NIST/TM1/Imager/';
savedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/NIST/TM'...
                    ,num2str(band),'/framedat/');                
%% H band

time1={'11-03-27';'11-04-44';'11-05-04';'11-05-27';'11-05-54'};
time2={'11-08-36';'11-08-59';'11-09-18';'11-09-35';'11-09-56'};
time3={'11-13-30';'11-13-51';'11-14-10';'11-14-44';'11-15-06'};
time4={'11-16-35';'11-16-59';'11-17-19';'11-17-40';'11-18-03'};
time_arr={time1,time2,time3,time4};

loaddir='/Volumes/HD1TB/CIBER/data/Cal_NIST/TM2/Imager/';
savedir=strcat('/Users/ytcheng/Documents/MATLAB/20160829_CalNIST/TM'...
                    ,num2str(band),'/framedat/');                
%%

for iset=1:numel(time_arr)
for itime=1:numel(time_arr{iset})
time=char(time_arr{iset}(itime));
files = dir(strcat(loaddir,'rawdat/*',time,'*')); 

frames=zeros(numel(files),1024,1024);
for ifr=1:numel(files)
infile = strcat(loaddir,'rawdat/',files(ifr).name);
frame = imrotate(fitsread(infile),270);    
frames(ifr,:,:)=frame; 
end

start_arr=zeros(1024);end_arr=zeros(1024);
for rr=1:1024
    for cc=1:1024
        ts=squeeze(frames(:,rr,cc));
        [fitstart,fitend]=reset_finder(ts);
        start_arr(rr,cc)=fitstart;
        end_arr(rr,cc)=fitend;
    end
end

starthist=histc(start_arr(:),1:numel(files));
endhist=histc(end_arr(:),1:numel(files));
[histost,indexst]=sort(starthist,'descend');
[histoed,indexed]=sort(endhist,'descend');
fts=indexst(1);sts=indexst(2);ffs=histost(1);sfs=histost(2);
fte=indexed(1);ste=indexed(2);ffe=histoed(1);sfe=histoed(2);
fps=ffs.*100./1024.^2;sps=sfs.*100./1024.^2;
fpe=ffe.*100./1024.^2;spe=sfe.*100./1024.^2;
disp(sprintf(strcat('%s,nfr=%d,start=%d(%f %%),end=%d(%f %%),',...
        '[start=%d(%f %%),end=%d(%f %%)]'),time,numel(files),fts,fps,...
        fte,fpe,sts,sps,ste,spe));

linefit=fastlinefit_frin(frames(fts:fte,:,:),0);

framedat.time=time;
framedat.data=frames(fts:fte,:,:);
framedat.map=linefit;
framedat.nfr=numel(files);
framedat.start=fts;framedat.end=fte;
framedat.start_freq=ffs;framedat.start_perc=fps;
framedat.end_freq=ffe;framedat.end_perc=fpe;
framedat.start2=sts;framedat.end2=ste;
framedat.start_freq2=sfs;framedat.start_perc2=sps;
framedat.end_freq2=sfe;framedat.end_perc2=spe;

save(strcat(savedir,time,'_framedat'),'framedat');
end
end  

