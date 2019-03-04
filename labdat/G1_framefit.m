%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%G1 cal with FF lab data.
%This code find resets to determine useful frames,
%and do the line fit to get map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band=1;
loaddir=strcat('/Users/ytcheng/ciber/data/Cal_FF/TM',num2str(band),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/FF/TM'...
    ,num2str(band),'/framedat/'); 

scanfile=dir(strcat(loaddir,'*0001.fits'));
for i=1:numel(scanfile)
time=scanfile(i).name(15:22);
files = dir(strcat(loaddir,'*',time,'*')); 

frames=zeros(numel(files),1024,1024);
for ifr=1:numel(files)
infile = strcat(loaddir,files(ifr).name);
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
disp(sprintf(strcat('%s,nfr=%d,start=%d(%f %%),end=%d(%f %%),',...
 'bad_start=%d, bad_end=%d'),time,numel(files),indexst(1),...
 histost(1).*100./1024.^2,indexed(1),histoed(1).*100./1024.^2,...
 numel(find(start_arr~=indexst(1))),numel(find(end_arr~=indexed(1)))));

mask=ones(1024);
mask(find(start_arr~=indexst(1)))=0;mask(find(end_arr~=indexed(1)))=0;


framedat.time=time;
framedat.data=frames(indexst(1):indexed(1),:,:);
framedat.nfr=numel(files);
framedat.start=indexst(1);framedat.end=indexed(1);
framedat.mask=mask;

% if available frames >=5, get linefit. 
if indexed(1)-indexst(1)>=4
linefit=fastlinefit_frin(frames(indexst(1):indexst(1)+4,:,:),0);
framedat.map5=linefit;
linefit=fastlinefit_frin(frames(indexst(1):indexst(1)+3,:,:),0);
framedat.map4=linefit;
linefit=fastlinefit_frin(frames(indexst(1):indexst(1)+2,:,:),0);
framedat.map3=linefit;
linefit=fastlinefit_frin(frames(indexst(1):indexst(1)+1,:,:),0);
framedat.map2=linefit;
else 
framedat.map5=[];framedat.map4=[];framedat.map3=[];framedat.map2=[];
end

save(strcat(savedir,time,'_framedat'),'framedat');
end