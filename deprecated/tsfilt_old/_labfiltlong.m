%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Filter with the long time stream
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
pixscale=7;
inst=1;
darkdir=strcat('/Volumes/HD1TB/CIBER/data/Dark_Raw/',...
    num2str(flight),'/TM',num2str(inst),'/');
[time_arr,start_arr,end_arr]=get_dark_long(flight,inst);
maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');

savedir='/Volumes/HD1TB/CIBER/tsfilt/sinfiltamp/DarkLong/';
%% do the filtering and save filt frames and map
%{
for i=1:numel(time_arr)
    disp(sprintf('idark=%d/%d',i,numel(time_arr)));

    time=time_arr{i};
    startfr=start_arr(i);
    endfr=end_arr(i);
    nfr=endfr-startfr+1;

    scanfile=dir(strcat(darkdir,'*',time,'*'));
    rawfrlong=zeros(nfr,1024,1024);
    for j=1:nfr
        fname=strcat(darkdir,scanfile(j+startfr-1).name);
        frame = imrotate(fitsread(fname),270);
        rawfrlong(j,:,:)=frame;
    end
    % Filtering with the full time stream
    [rawmap,rawoff]=linfit_map(rawfrlong,'verbose',0); 
    [~,maskin]=get_skymap(rawmap,mask_inst,2,5);
    [filtfrlong] = imager_filtts(rawfrlong,'mapin',rawmap,...
     'maskin',maskin,'offin',rawoff,'verbose',0,'makeplot',0);
    for j=1:nfr
        frame=squeeze(filtfrlong(j,:,:));
        frind=j+startfr-1;
        save(sprintf('%s/frame/dark%s_%.3d',savedir,time,frind),'frame');
    end  
    filtmap=linfit_map(filtfrlong,'verbose',0);
    save(sprintf('%s/map/dark%s_%.3d_%.3d',savedir,time,startfr,endfr),...
                                                              'filtmap');
end
%}
%% stack DC template 
dcstack=zeros(1024);
maskstack=zeros(1024);

for i=1:numel(time_arr)
    time=time_arr{i};
    startfr=start_arr(i);
    endfr=end_arr(i);
    nfr=endfr-startfr+1;
    wi=nfr*(nfr^2-1);
    load(sprintf('%s/map/dark%s_%.3d_%.3d',savedir,time,startfr,endfr),...
                                                              'filtmap');
    [~,maski]=get_skymap(filtmap,mask_inst,4,5);
    dcstack=dcstack+(filtmap.*maski.*wi);
    maskstack=maskstack+maski.*wi;
end
dcstack=dcstack./maskstack;dcstack(find(dcstack~=dcstack))=0;
dcstack(find(dcstack==inf))=0;dcstack(find(dcstack==-inf))=0;
DCtemplate=dcstack;
save(sprintf('%s/DCtemplate',savedir),'DCtemplate');