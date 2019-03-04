date='05-13-2013';
band=2;

datadir=strcat('/Users/ytcheng/ciber/data/focus/',date,'/');
savedatdir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/',...
            'Focus/slopedata/',date,'/TM',num2str(band),'/'); 
saveplotdir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/',...
            'Focus/slopemap/',date,'/TM',num2str(band),'/'); 
% only select data with nfr>=2, so reset_finder can goes on
scantime=dir(strcat(datadir,'C',num2str(band),'*0002.fits'));
%%
for itime=1:numel(scantime)
time=scantime(itime).name(15:22);
scanframe=dir(strcat(datadir,scantime(itime).name(1:end-9),'*'));

% write frames arr
frames=zeros(numel(scanframe),1024,1024);
for ifr=1:numel(scanframe)
    fname=strcat(datadir,scanframe(ifr).name);
    frame = imrotate(fitsread(fname),270);
    frames(ifr,:,:)=frame;
end

% find resets
start_arr=zeros(1024);end_arr=zeros(1024);
for rr=1:1024
    for cc=1:1024
        ts=squeeze(frames(:,rr,cc));
        [fitstart,fitend]=reset_finder(ts);
        start_arr(rr,cc)=fitstart;
        end_arr(rr,cc)=fitend;
    end
end
% reset stat
starthist=histc(start_arr(:),1:numel(scanframe));
endhist=histc(end_arr(:),1:numel(scanframe));
[histost,indexst]=sort(starthist,'descend');
[histoed,indexed]=sort(endhist,'descend');
disp(sprintf(strcat('%d,%s,nfr=%d,start=%d(%f %%),end=%d(%f %%),',...
 'bad_start=%d, bad_end=%d'),itime,time,numel(scanframe),indexst(1),...
 histost(1).*100./1024.^2,indexed(1),histoed(1).*100./1024.^2,...
 numel(find(start_arr~=indexst(1))),numel(find(end_arr~=indexed(1)))));

% mask bad reset stat pix
mask=ones(1024);
mask(find(start_arr~=indexst(1)))=0;mask(find(end_arr~=indexed(1)))=0;

if indexed(1)-indexst(1)+1>2
    framedat.band=band;
    framedat.date=date;
    framedat.time=time;
    framedat.nfr=numel(scanframe);
    framedat.data=frames(indexst(1):indexed(1),:,:);
    framedat.start=indexst(1);framedat.end=indexed(1);
    framedat.mask=mask;

    framedat.nfr_arr=2:indexed(1)-indexst(1)+1;
    map_arr=zeros(numel(framedat.nfr_arr),1024,1024);
    for infr=2:indexed(1)-indexst(1)+1
        map=linfit_map(frames(indexst(1):indexst(1)+infr-1,:,:),...
            'verbose',0);map_arr(infr-1,:,:)=map;
    end
    framedat.map_arr=map_arr;
    save(strcat(savedatdir,time,'_framedat'),'framedat');
    %plot the longest int time map 
    imageclip(map);
    plotname=strcat(saveplotdir,scantime(itime).name(1:end-10));
    print(plotname,'-dpng');close
else
    disp(sprintf('%s,NOT WRITE',time));
end

end
