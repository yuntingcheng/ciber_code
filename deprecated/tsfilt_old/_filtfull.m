%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filter lab dark and flgiht
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
pixscale=7;

darkdir=strcat('/Volumes/HD1TB/CIBER/data/Dark_Raw/40030/TM',...
    num2str(inst),'/');
maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');
%%
for ifield=8:-1:1
savedir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
bigmask=alldat(ifield).bigmask.*mask_inst;
clear alldat
dt=get_dark_times(flight,inst,ifield);
nfrhalf=dt.nfrhalf;
nfrfull=dt.nfr;

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
     'offin',rawoff,'verbose',0,'makeplot',0);
    
    labdat(i).filtfr=filtfr;    
end

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
    
    labdat(i).rawmap1=rawmap1;
    labdat(i).rawmap2=rawmap2;
    labdat(i).rawmapf=rawmapf;
    labdat(i).filtmap1=filtmap1;
    labdat(i).filtmap2=filtmap2;
    labdat(i).filtmapf=filtmapf;

end

labdat1=rmfield(labdat,{'rawfr','rawmap','rawoff','filtfr'});
clear labdat

for i=1:numel(dt.time)
    labdat.rawmap1=labdat1(i).rawmap1;
    labdat.rawmap2=labdat1(i).rawmap2;
    labdat.rawmapf=labdat1(i).rawmapf;
    labdat.filtmap1=labdat1(i).filtmap1;
    labdat.filtmap2=labdat1(i).filtmap2;
    labdat.filtmapf=labdat1(i).filtmapf;

    save(strcat(savedir,'labdat',num2str(i)),'labdat');
end

%%% filter flight
fields=get_fields(flight,inst);
name=fields(ifield).name;
[rawfr] = get_data_frames (inst,name,'flight',flight,'verbose',0);
rawfr=rawfr(3:end,:,:);

%%% line fit rawmap
disp(sprintf('linfit rawmap'));
[rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 
[~,maskin1]=get_skymap(rawmap,maskin,2,5);

%%% do filtering
disp(sprintf('start filtering'));
[filtfr] = imager_filtts(rawfr,'mapin',rawmap,'maskin',maskin1,...
    'offin',rawoff,'verbose',0,'makeplot',0);

%%% flight diff maps %%%%%

disp(sprintf('linefit flight diff'));

fr1=rawfr(1:nfrhalf,:,:);
fr2=rawfr(nfrhalf+1:2*nfrhalf,:,:);
[rawmap1] =linfit_map(fr1,'verbose',0);
[rawmap2] =linfit_map(fr2,'verbose',0);
[rawmapf] =linfit_map(rawfr,'verbose',0); 

fr1=filtfr(1:nfrhalf,:,:);
fr2=filtfr(nfrhalf+1:2*nfrhalf,:,:);
[filtmap1] =linfit_map(fr1,'verbose',0);
[filtmap2] =linfit_map(fr2,'verbose',0);
[filtmapf] =linfit_map(filtfr,'verbose',0);
flightmap.rawmap1=rawmap1;
flightmap.rawmap2=rawmap2;
flightmap.rawmapf=rawmapf;
flightmap.filtmap1=filtmap1;
flightmap.filtmap2=filtmap2;
flightmap.filtmapf=filtmapf;
save(strcat(savedir,'flightmap'),'flightmap');
end

