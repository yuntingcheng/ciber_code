%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Do the line fit on dark current time stream data to get DC map.
% - Do both full frames, and half frames.
% - The number of frames for each field, the good frame
% range DarkList.xlsx, and coded in get_dark_times.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=36265;
band=1;

switch flight
case {36265,36277}
darkdir=strcat('/Users/ytcheng/ciber/data/Dark_Raw/362xx/TM',...
    num2str(band),'/');
case 40030
darkdir=strcat('/Users/ytcheng/ciber/data/Dark_Raw/40030/TM',...
    num2str(band),'/');
end

savedir=strcat('/Users/ytcheng/ciber/doc/20150810_DarkCurrent/',...
    num2str(flight),'/');
%%
for field=1:numel(get_fields(flight,band))
dt=get_dark_times(flight,band,field);

%%%creat folders%%%
foldname=sprintf('%sTM%d/%s',savedir,band,dt.name);
if exist(foldname,'dir');rmdir(foldname,'s');end
if ~exist(foldname, 'dir');mkdir(foldname);end

foldname=sprintf('%sTM%d/%s/full/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end
foldname=sprintf('%sTM%d/%s/full/data/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end
foldname=sprintf('%sTM%d/%s/full/plot/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end

foldname=sprintf('%sTM%d/%s/first/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end
foldname=sprintf('%sTM%d/%s/first/data/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end
foldname=sprintf('%sTM%d/%s/first/plot/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end

foldname=sprintf('%sTM%d/%s/second/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end
foldname=sprintf('%sTM%d/%s/second/data/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end
foldname=sprintf('%sTM%d/%s/second/plot/',savedir,band,dt.name);
if ~exist(foldname, 'dir');mkdir(foldname);end

%%%%%%%%%%%%%%%%% full %%%%%%%%%%%%%%%%%%%%%%%%
%%% do line fit %%%
for i=1:numel(dt.time)
time=dt.time{i};
scanfile=dir(strcat(darkdir,'*',time,'*'));
frames=zeros(dt.nfr,1024,1024);
for j=dt.frdown(i):dt.frup(i)
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    frames(j-dt.frdown(i)+1,:,:)=frame;
end
darkmap=linfit_map(frames,'verbose',0);
%%% save data %%%
dataname=sprintf('%sTM%d/%s/full/data/dark_%s_%d_%d',...
                    savedir,band,dt.name,time,dt.frdown(i),dt.frup(i));
save(dataname,'darkmap');
%%% save plot %%%
imageclip(darkmap);
plotname=sprintf('%sTM%d/%s/full/plot/dark_%s_%d_%d',...
                    savedir,band,dt.name,time,dt.frdown(i),dt.frup(i));
title(sprintf('TM%d %s %d-%d',band,time,dt.frdown(i),dt.frup(i)))
print(plotname,'-dpng');close

pr=sprintf('%s,nfr=%d,full_%s,%d-%d',...
                dt.name,dt.nfr,time,dt.frdown(i),dt.frup(i));disp(pr);
end
%%%%%%%%%%%%%%%%% first %%%%%%%%%%%%%%%%%%%%%%%%
%%% do line fit %%%
for i=1:numel(dt.time)
time=dt.time{i};
scanfile=dir(strcat(darkdir,'*',time,'*'));
frames=zeros(dt.nfrhalf,1024,1024);
for j=dt.frdown1(i):dt.frup1(i)
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    frames(j-dt.frdown1(i)+1,:,:)=frame;
end
darkmap=linfit_map(frames,'verbose',0);
%%% save data %%%
dataname=sprintf('%sTM%d/%s/first/data/dark_%s_%d_%d',...
                    savedir,band,dt.name,time,dt.frdown1(i),dt.frup1(i));
save(dataname,'darkmap');
%%% save plot %%%
imageclip(darkmap);
plotname=sprintf('%sTM%d/%s/first/plot/dark_%s_%d_%d',...
                    savedir,band,dt.name,time,dt.frdown1(i),dt.frup1(i));
title(sprintf('TM%d %s %d-%d',band,time,dt.frdown1(i),dt.frup1(i)))
print(plotname,'-dpng');close

pr=sprintf('%s,nfrhalf=%d,first_%s,%d-%d',...
              dt.name,dt.nfrhalf,time,dt.frdown1(i),dt.frup1(i));disp(pr);
end
%%%%%%%%%%%%%%%%% second %%%%%%%%%%%%%%%%%%%%%%%%
%%% do line fit %%%
for i=1:numel(dt.time)
time=dt.time{i};
scanfile=dir(strcat(darkdir,'*',time,'*'));
frames=zeros(dt.nfrhalf,1024,1024);
for j=dt.frdown2(i):dt.frup2(i)
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    frames(j-dt.frdown2(i)+1,:,:)=frame;
end
darkmap=linfit_map(frames,'verbose',0);
%%% save data %%%
dataname=sprintf('%sTM%d/%s/second/data/dark_%s_%d_%d',...
                    savedir,band,dt.name,time,dt.frdown2(i),dt.frup2(i));
save(dataname,'darkmap');
%%% save plot %%%
imageclip(darkmap);
plotname=sprintf('%sTM%d/%s/second/plot/dark_%s_%d_%d',...
                    savedir,band,dt.name,time,dt.frdown2(i),dt.frup2(i));
title(sprintf('TM%d %s %d-%d',band,time,dt.frdown2(i),dt.frup2(i)))
print(plotname,'-dpng');close

pr=sprintf('%s,nfrhalf=%d,second_%s,%d-%d',...
              dt.name,dt.nfrhalf,time,dt.frdown2(i),dt.frup2(i));disp(pr);
end

end


