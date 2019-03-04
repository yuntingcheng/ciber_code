%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Do the line fit on dark current time stream data to get DC map.
% - Do both full frames, and half frames.
% - The number of frames for each field, the good frame
% range DarkList.xlsx, and coded in get_dark_times.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
band=1;

darkdir=strcat('/Volumes/HD1TB/CIBER/data/Dark_Raw/40030/TM',...
    num2str(band),'/');

savedir=strcat('/Users/ytcheng/ciber/doc/20150810_DarkCurrent/',...
    '40030int/');
%% dark data
dt=get_dark_times(flight,band,8);
for infr=2:13
    darkmap1_arr=zeros(numel(dt.time),1024,1024);
    darkmap2_arr=zeros(numel(dt.time),1024,1024);
    
    for i=1:numel(dt.time)
        time=dt.time{i};
        scanfile=dir(strcat(darkdir,'*',time,'*'));

        frames=zeros(infr,1024,1024);
        for j=dt.frdown1(i):dt.frdown1(i)+infr-1
            fname=strcat(darkdir,scanfile(j).name);
            frame = imrotate(fitsread(fname),270);
            frames(j-dt.frdown1(i)+1,:,:)=frame;
        end
        darkmap1=linfit_map(frames,'verbose',0);
        

        frames=zeros(infr,1024,1024);
        for j=dt.frdown2(i):dt.frdown2(i)+infr-1
            fname=strcat(darkdir,scanfile(j).name);
            frame = imrotate(fitsread(fname),270);
            frames(j-dt.frdown2(i)+1,:,:)=frame;
        end
        darkmap2=linfit_map(frames,'verbose',0);
        
        darkmap1_arr(i,:,:)=darkmap1;
        darkmap2_arr(i,:,:)=darkmap2;
        
        disp(sprintf('nfr=%d,dark %d/%d',infr,i,numel(dt.time)));
    end
    rnintdat(infr).darkmap1=darkmap1_arr;
    rnintdat(infr).darkmap2=darkmap2_arr;
end

save(sprintf('%s/TM%d_rnintdat',savedir,band),'rnintdat');
%% flight data

ft = get_field_times(flight,band);
fields=get_fields(flight,band);
for i=4:8
name=fields(i).name;
[frames] = get_data_frames(band,name,'flight',flight,'verbose',0);

frdat(i).frames=frames;    
end

for infr=2:13
    flightmap1_arr=zeros(8,1024,1024);
    flightmap2_arr=zeros(8,1024,1024);
    
    for i=4:8
    frames=frdat(i).frames;

    frames=frames(3:end,:,:);
    nfrhalf=floor(size(frames,1)/2);

    [rawmap1]=linfit_map(frames(1:infr,:,:),'verbose',0);
    [rawmap2]=linfit_map(frames(1+nfrhalf:nfrhalf+infr,:,:),'verbose',0);
    flightmap1_arr(i,:,:)=rawmap1;
    flightmap2_arr(i,:,:)=rawmap2;
    disp(sprintf('nfr=%d,field %d',infr,i));
    end
    flintdat(infr).flightmap1=flightmap1_arr;
    flintdat(infr).flightmap1=flightmap2_arr;
end
save(sprintf('%s/TM%d_flintdat',savedir,band),'flintdat');
