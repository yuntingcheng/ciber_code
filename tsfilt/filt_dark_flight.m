function filt_dark_flight(flight,inst,ifield)
if inst==1
    sin_freq=9.503;
elseif inst==2
    sin_freq=9.538;
end
mypaths=get_paths(flight);
darkdir=strcat(mypaths.Dark_Raw,'TM',num2str(inst),'/');
maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');
dt=get_dark_times(flight,inst,ifield);
nfrhalf=dt.nfrhalf;
nfr_arr=2:dt.nfrhalf;

% make dir for saving the result
savedir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
if ~exist(savedir, 'dir')
    mkdir(savedir);
end

% get mask from MZ dr
load(sprintf('%sTM%d_%s_dr150206.mat',mypaths.release,inst,dt.name));
maskin=double(~data.mask.mask);
maskin=maskin.*mask_inst;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% filter the lab dark %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get raw frames, rawmap, rawoff, and combined mask
for i=1:numel(dt.time)
    
    % read in raw frames   
    disp(sprintf('retrive frames for ifield=%d,idark=%d/%d',...
        ifield,i,numel(dt.time)));
    scanfile=dir(strcat(darkdir,'*',dt.time{i},'*'));
    rawfr=zeros(dt.nfr,1024,1024);
    for j=dt.frdown(i):dt.frup(i)
        fname=strcat(darkdir,scanfile(j).name);
        frame = imrotate(fitsread(fname),270);
        rawfr(j-dt.frdown(i)+1,:,:)=frame;
    end

    %%% line fit rawmap
    [rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 
    
    %%% merge mask
    [~,maskin]=get_skymap(rawmap,maskin,4,5);
    
    
    labdat(i).rawfr=rawfr;
    labdat(i).rawmap=rawmap;
    labdat(i).rawoff=rawoff;
end

save(strcat(savedir,'maskin'),'maskin');

%%% filter lab dark
for i=1:numel(dt.time)
    disp(sprintf('filter ifield=%d,idark=%d/%d',...
        ifield,i,numel(dt.time)));
    
    rawfr=labdat(i).rawfr;
    rawmap=labdat(i).rawmap;
    rawoff=labdat(i).rawoff;
    
    [~,maskin1]=get_skymap(rawmap,maskin,2,5);
    
    %%% do filtering
    %[filtfr] = imager_filtts(rawfr,rawoff,maskin1,'inst',inst,...
    % 'sin_freq',sin_freq,'verbose',1,'makeplot',1);
    [filtfr] = imager_filtts(rawfr,rawoff,maskin1,'inst',inst,...
    'donotch',1,'dosinfilt',0,'verbose',0,'makeplot',0,...
    'notch_nus',[sin_freq],'notch_width',[0.1],'notch_iter',1); 

    labdat(i).filtfr=filtfr;  
    
end
%%
%%% linefit diff map
for i=1:numel(dt.time)
    disp(sprintf('linefit,idark=%d/%d',i,numel(dt.time)));
   
    fr1=labdat(i).rawfr(1:nfrhalf,:,:);
    fr2=labdat(i).rawfr(nfrhalf+1:2*nfrhalf,:,:);        
    [rawmap1] =linfit_map(fr1,'verbose',0);
    [rawmap2] =linfit_map(fr2,'verbose',0);
    [rawmapf] =linfit_map(labdat(i).rawfr,'verbose',0); 

    fr1=labdat(i).filtfr(1:nfrhalf,:,:);
    fr2=labdat(i).filtfr(nfrhalf+1:2*nfrhalf,:,:);        
    [filtmap1] =linfit_map(fr1,'verbose',0);
    [filtmap2] =linfit_map(fr2,'verbose',0);
    [filtmapf] =linfit_map(labdat(i).filtfr,'verbose',0); 
    
    labmap.rawmap1=rawmap1;
    labmap.rawmap2=rawmap2;
    labmap.rawmapf=rawmapf;
    labmap.filtmap1=filtmap1;
    labmap.filtmap2=filtmap2;
    labmap.filtmapf=filtmapf;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rawmap_arr=zeros(numel(nfr_arr),1024,1024);
    filtmap_arr=zeros(numel(nfr_arr),1024,1024);
    for infr=1:numel(nfr_arr)
        nfr=nfr_arr(infr);
        disp(sprintf('linefit,idark=%d/%d,nfr=%d',i,numel(dt.time),nfr));
   
        fr1=labdat(i).rawfr(1:nfr,:,:);
        fr2=labdat(i).rawfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);        
        [map1] =linfit_map(fr1,'verbose',0);
        [map2] =linfit_map(fr2,'verbose',0);        
        diffraw=(map1-map2)./2;
        rawmap_arr(infr,:,:)=diffraw;
        
        fr1=labdat(i).filtfr(1:nfr,:,:);
        fr2=labdat(i).filtfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);        
        [map1] =linfit_map(fr1,'verbose',0);
        [map2] =linfit_map(fr2,'verbose',0);        
        difffilt=(map1-map2)./2;
        filtmap_arr(infr,:,:)=difffilt;        
    end
    labmap.rawmap_arr=rawmap_arr;
    labmap.filtmap_arr=filtmap_arr;
    
    save(strcat(savedir,'labmap',num2str(i)),'labmap');
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% filter the flight %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rawfr] = get_data_frames(inst,dt.name,'flight',flight,'verbose',0);
rawfr=rawfr(3:end,:,:);

disp(sprintf('linfit flight rawmap'));
[rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 
[~,maskin1]=get_skymap(rawmap,maskin,2,5);

disp(sprintf('start filtering flight'));
%[filtfr] = imager_filtts(rawfr,rawoff,maskin1,'inst',inst,...
%    'sin_freq',sin_freq,'verbose',0,'makeplot',0);
[filtfr] = imager_filtts(rawfr,rawoff,maskin1,'inst',inst,...
    'donotch',1,'dosinfilt',0,'verbose',0,'makeplot',0,...
    'notch_nus',[sin_freq],'notch_width',[0.1],'notch_iter',1); 

%flight map line fit
fr1=rawfr(1:nfrhalf,:,:);
fr2=rawfr(nfrhalf+1:2*nfrhalf,:,:);
[rawmap1] =linfit_map(fr1,'verbose',0);
[rawmap2] =linfit_map(fr2,'verbose',0);
[rawmapf] =linfit_map(rawfr,'verbose',0); 
flightmap.rawmap1=rawmap1;
flightmap.rawmap2=rawmap2;
flightmap.rawmapf=rawmapf;

fr1=filtfr(1:nfrhalf,:,:);
fr2=filtfr(nfrhalf+1:2*nfrhalf,:,:);
[filtmap1] =linfit_map(fr1,'verbose',0);
[filtmap2] =linfit_map(fr2,'verbose',0);
[filtmapf] =linfit_map(filtfr,'verbose',0);
flightmap.filtmap1=filtmap1;
flightmap.filtmap2=filtmap2;
flightmap.filtmapf=filtmapf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawmap_arr=zeros(numel(nfr_arr),1024,1024);
filtmap_arr=zeros(numel(nfr_arr),1024,1024);

for infr=1:numel(nfr_arr)
    nfr=nfr_arr(infr);
    disp(sprintf('linefit flight diff,nfr=%d',nfr));

    fr1=rawfr(1:nfr,:,:);
    fr2=rawfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);
    [map1] =linfit_map(fr1,'verbose',0);
    [map2] =linfit_map(fr2,'verbose',0);
    diff=(map1-map2)./2;
    rawmap_arr(infr,:,:)=diff;
    
    fr1=filtfr(1:nfr,:,:);
    fr2=filtfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);
    [map1] =linfit_map(fr1,'verbose',0);
    [map2] =linfit_map(fr2,'verbose',0);
    diff=(map1-map2)./2;
    filtmap_arr(infr,:,:)=diff;
end
flightmap.nfr_arr=nfr_arr;
flightmap.rawmap_arr=rawmap_arr;
flightmap.filtmap_arr=filtmap_arr;

save(strcat(savedir,'flightmap'),'flightmap');

return 