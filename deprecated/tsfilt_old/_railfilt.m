%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Do the line fit on dark current time stream data to get DC map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
pixscale=7;
fields=get_fields(flight,inst);

maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');
%%
ifield=8;
%savedir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
%            'sinfilt/DiffMap/field',num2str(ifield),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170209_TsFilter/',...
        'sinfiltrange/DiffMap/field',num2str(ifield),'/');

load(strcat(savedir,'maskin'),'maskin');
dt=get_dark_times(flight,inst,ifield);
%% filter rail dark
load('/Users/ytcheng/ciber/data/raildark_40030',...
strcat('fr',num2str(inst)));%37x1024x1024 data
rawfr=fr1;
nfr=size(rawfr,1);
nfrhalf=floor(nfr/2);

%%% line fit rawmap
disp(sprintf('linfit rawmap'));
[rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 

%[~,maskin]=get_skymap(rawmap,maskin,4,5);
[~,maskin]=get_skymap(rawmap,maskin,2,5);

railmask=maskin;
%save(strcat(savedir,'railmask'),'railmask');

%%% do filtering
disp(sprintf('start filtering'));
[filtfr] = imager_filtts(rawfr,'mapin',rawmap,...
'maskin',maskin,'offin',rawoff,'verbose',1,'makeplot',1);
    
%%% line fit filtmap
disp(sprintf('linfit filtmap'));
[filtmap] =linfit_map(filtfr,'verbose',0);

raildat.rawfr=rawfr;
raildat.filtfr=filtfr;
raildat.maskin=maskin;
raildat.nfr=nfr;
raildat.nfrhalf=nfrhalf;
raildat.rawmap=rawmap;
raildat.filtmap=filtmap;
%%% do diff for rail (raw & flight)
nfr_arr=2:dt.nfrhalf;

rawmap_arr=zeros(numel(nfr_arr),1024,1024);
filtmap_arr=zeros(numel(nfr_arr),1024,1024);

for infr=1:numel(nfr_arr)
    nfr=nfr_arr(infr);
    disp(sprintf('linefit rail diff,nfr=%d',nfr));

    fr1=raildat.rawfr(1:nfr,:,:);
    fr2=raildat.rawfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);
    [map1] =linfit_map(fr1,'verbose',0);
    [map2] =linfit_map(fr2,'verbose',0);
    diff=(map1-map2)./2;
    rawmap_arr(infr,:,:)=diff;
    
    fr1=raildat.filtfr(1:nfr,:,:);
    fr2=raildat.filtfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);
    [map1] =linfit_map(fr1,'verbose',0);
    [map2] =linfit_map(fr2,'verbose',0);
    diff=(map1-map2)./2;
    filtmap_arr(infr,:,:)=diff;
end
railmap.rawmap_arr=rawmap_arr;
railmap.filtmap_arr=filtmap_arr;
save(strcat(savedir,'railmap'),'railmap');
