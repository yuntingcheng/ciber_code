function fitdat=fit_focus_g1(framedir,mask_inst,nfr,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process the focus map with slope data in 'framedir'.
% Only use the ambient, so need a man_mask out the source. The source mask
% is generate by smoothing the map, and mask the brightest mask_perc
% percent of the smoothed map
% input:
% - framedir: where the input slope data is 
% - mask_inst: instrument mask
% - nfr: # of frames to extract from slope data
% - sm_sig(Optional): sigma used in smooth the map (to get man_mask)
% - mask_perc(Optional): percentage of pix to be masked in smoothed map
% - quad_corr(Optional):(1/0)mask top right quadrant or not
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse data
  p = inputParser;
  p.addRequired('framedir',@ischar);
  p.addRequired('mask_inst',@isnumeric);
  p.addRequired('nfr',@isnumeric);
  p.addOptional('sm_sig',5,@isnumeric);
  p.addOptional('mask_perc',30,@isnumeric);
  p.addOptional('quad_corr',0,@isnumeric);
  p.parse(framedir,mask_inst,nfr,varargin{:});
  framedir     = p.Results.framedir;
  mask_inst  = p.Results.mask_inst;
  nfr  = p.Results.nfr;
  sm_sig  = p.Results.sm_sig;
  mask_perc  = p.Results.mask_perc;
  quad_corr  = p.Results.quad_corr;
  clear p varargin;
  

%%% scan the dir  
scanfile=dir(strcat(framedir,'*.mat'));

% count # of map with desired nfr
Nmap=0;
for i=1:numel(scanfile)    
    load(strcat(framedir,scanfile(i).name),'framedat');
    if  find(framedat.nfr_arr==nfr);Nmap=Nmap+1;end
end

% if Nmap less than 2=> skip the rest process
if Nmap<2
disp(sprintf('%s, no nfr=%d data to process',framedir,nfr));fitdat=[];

% keep going if Nmap>=2
else

% write maps and sigma-clip mask with given nfr into map_arr & mask_arr
map_arr=zeros(Nmap,1024,1024);
mask_arr=zeros(Nmap,1024,1024);
count=0;
for i=1:numel(scanfile) 
    load(strcat(framedir,scanfile(i).name),'framedat');
    if find(framedat.nfr_arr==nfr)
        count=count+1;
        map=squeeze(framedat.map_arr(find(framedat.nfr_arr==nfr),:,:));
        map_arr(count,:,:)=map;
        
        mask_quad=ones(1024);
        if quad_corr; mask_quad(513:1024,513:1024)=0;end
        
        [~,mask]=get_skymap(map,mask_inst.*mask_quad,3);
        mask_arr(count,:,:)=mask;
    end
end

% get man_mask
mapstack=zeros(1024);maskstack=zeros(1024);
for i=1:Nmap
    map=squeeze(map_arr(i,:,:));
    mask=squeeze(mask_arr(i,:,:));    
    mapstack=mapstack+map.*mask;
    maskstack=maskstack+mask;
end
mapstack=mapstack./maskstack;mapstack(find(mapstack~=mapstack))=0;
maskstack(find(maskstack))=1;
sm_map=fillpadsmooth(mapstack,maskstack,sm_sig);
a=sm_map(find(maskstack));
man_mask=ones(1024);man_mask(find(sm_map<prctile(a,mask_perc)))=0;

% exclude map with crazy masked mean (mean > 2sigma)
mean_arr=zeros(1,Nmap);
for i=1:Nmap    
    map=squeeze(map_arr(i,:,:));
    mask=squeeze(mask_arr(i,:,:));    
    mm=map.*mask.*man_mask;
    mean_arr(i)=mean(mm(find(mm)));    
end
bad_ind=find(mean_arr>median(mean_arr)+std(mean_arr)*2 | ...
     mean_arr<median(mean_arr)-std(mean_arr)*2 );
use_ind=1:Nmap;
use_ind(bad_ind)=[];
mean_arr(bad_ind)=[];
[~,sortorder]=sort(mean_arr);
use_ind=use_ind(sortorder);

% throw out last one if we have odd number of data
if mod(numel(use_ind),2)==1;use_ind(end)=[];end

% calculate the mean and var of the map
avg_arr=zeros(1,numel(use_ind)/2);
var_arr=zeros(1,numel(use_ind)/2);
diffmap_arr=zeros(numel(use_ind)/2,1024,1024);

for i=1:numel(use_ind)/2
    map1=squeeze(map_arr(use_ind(i),:,:));
    mask1=squeeze(mask_arr(use_ind(i),:,:));    
    mm1=map1.*mask1.*man_mask;
    mean1=mean(mm1(find(mm1)));
    
    map2=squeeze(map_arr(use_ind(i+1),:,:));
    mask2=squeeze(mask_arr(use_ind(i+1),:,:));    
    mm2=map2.*mask2.*man_mask;
    mean2=mean(mm2(find(mm2)));
    
    diff=(mm1-mm2).*mask1.*mask2.*man_mask./sqrt(2);  
    
    avg_arr(i)=(mean1+mean2)/2;
    var_arr(i)=var(diff(find(diff)));
    diffmap_arr(i,:,:)=diff;
end

% find bad map: exclude var > min(var)*1.05 (more than 5%)
badindex=find(var_arr>min(var_arr).*1.05);
avg_arr(badindex)=[];var_arr(badindex)=[];diffmap_arr(badindex,:,:)=[];

% save results
fitdat.framedir=framedir;
fitdat.nfr=nfr;
fitdat.sm_sig=sm_sig;
fitdat.mask_perc=mask_perc;
fitdat.diffmap_arr=diffmap_arr;
fitdat.man_mask=man_mask;
fitdat.avg_arr=avg_arr;
fitdat.var_arr=var_arr;

end
return
