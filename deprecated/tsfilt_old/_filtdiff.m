%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter lab dark and flgiht
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;

darkdir=strcat('/Volumes/HD1TB/CIBER/data/Dark_Raw/40030/TM',...
    num2str(inst),'/');
maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');
%%
for ifield=1:8
savedir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/DiffMap/field',num2str(ifield),'/');
plotdir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
       'sinfiltamp/plot/field',num2str(ifield),'/');

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
bigmask=alldat(ifield).bigmask.*mask_inst;
clear alldat
dt=get_dark_times(flight,inst,ifield);
nfrhalf=dt.nfrhalf;
%%% get raw frames, rawmap, rawoff, and combined mask
maskin=bigmask;
for i=1:numel(dt.time)
    
    %%% read in raw frames   
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
    disp(sprintf('linfit rawmap'));
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
    
    [~,maskin1]=get_skymap(rawmap,bigmask,2,5);
    
    %%% do filtering
    [filtfr] = imager_filtts(rawfr,'mapin',rawmap,'maskin',maskin1,...
     'offin',rawoff,'verbose',0,'makeplot',0,'plotdir',plotdir);
    
    labdat(i).filtfr=filtfr;    
end

%%% linefit diff map
nfr_arr=2:dt.nfrhalf;
for i=1:numel(dt.time)
    rawmap_arr=zeros(numel(nfr_arr),1024,1024);
    filtmap_arr=zeros(numel(nfr_arr),1024,1024);
    for infr=1:numel(nfr_arr)
        nfr=nfr_arr(infr);
        disp(sprintf('linefit,idark=%d/%d,nfr=%d',...
        i,numel(dt.time),nfr));
   
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
    labdat(i).rawmap_arr=rawmap_arr;
    labdat(i).filtmap_arr=filtmap_arr;

end

%%%
labdat1=rmfield(labdat,{'rawfr','rawmap','rawoff','filtfr'});
clear labdat
for i=1:numel(dt.time)
    labdat.rawmap_arr=labdat1(i).rawmap_arr;
    labdat.filtmap_arr=labdat1(i).filtmap_arr;
    save(strcat(savedir,'labdat',num2str(i)),'labdat');
end



%%% filter flight
fields=get_fields(flight,inst);
name=fields(ifield).name;
[rawfr] = get_data_frames (inst,name,'flight',flight,'verbose',0);
rawfr=rawfr(3:end,:,:);
nfr=size(rawfr,1);
nfrhalf=floor(nfr/2);

%%% line fit rawmap
disp(sprintf('linfit rawmap'));
[rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 
[~,maskin1]=get_skymap(rawmap,maskin,2,5);

%%% do filtering
disp(sprintf('start filtering'));
[filtfr] = imager_filtts(rawfr,'mapin',rawmap,'maskin',maskin1,...
    'offin',rawoff,'verbose',0,'makeplot',0,'plotdir',plotdir);

%%% flight diff maps %%%%%
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

flightmap.rawmap_arr=rawmap_arr;
flightmap.filtmap_arr=filtmap_arr;
save(strcat(savedir,'flightmap'),'flightmap');
end