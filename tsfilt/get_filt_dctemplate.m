function get_filt_dctemplate(flight,inst)
if inst==1
    sin_freq=9.503;
elseif inst==2
    sin_freq=9.538;
end
mypaths=get_paths(flight);
darkdir=strcat(mypaths.Dark_Raw,'TM',num2str(inst),'/');

[time_arr,start_arr,end_arr]=get_dark_long(flight,inst);
maskinstdir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
load(strcat(maskinstdir,'band',num2str(inst),'_mask_inst'),'mask_inst');
savedir=sprintf('%sTM%d/',mypaths.filtmap,inst);
%% do the filtering and save filt frames and map

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
    %[filtfrlong] = imager_filtts(rawfrlong,'mapin',rawmap,'maskin',...
    % maskin,'offin',rawoff,'sin_freq',sin_freq,'verbose',0,'makeplot',0);
    [filtfrlong] = imager_filtts(rawfrlong,rawoff,maskin,'inst',inst,...
    'donotch',1,'dosinfilt',0,'verbose',0,'makeplot',0,...
    'notch_nus',[sin_freq],'notch_width',[0.1],'notch_iter',1); 

    filtmap=linfit_map(filtfrlong,'verbose',0);
    
    % save the result
    darklongdat.mask_inst=mask_inst;
    darklongdat.map(i).rawmap=rawmap;
    darklongdat.map(i).filtmap=filtmap;
    darklongdat.map(i).nfr=nfr;
end

%% stack DC template 
dcstack=zeros(1024);
maskstack=zeros(1024);

for i=1:numel(time_arr)
    nfr=darklongdat.map(i).nfr;
    filtmap=darklongdat.map(i).filtmap;
    wi=nfr*(nfr^2-1);
    [~,maski]=get_skymap(filtmap,darklongdat.mask_inst,4,5);
    [filtmap]=dc_offset_remove(filtmap,maski);
    dcstack=dcstack+(filtmap.*maski.*wi);
    maskstack=maskstack+maski.*wi;
end
dcstack=dcstack./maskstack;dcstack(find(dcstack~=dcstack))=0;
dcstack(find(dcstack==inf))=0;dcstack(find(dcstack==-inf))=0;
darklongdat.DCtemplate=dcstack;
%% stack DC template with no mask_inst, for darkstat generate
dcstack=zeros(1024);
maskstack=zeros(1024);

for i=1:numel(time_arr)
    nfr=darklongdat.map(i).nfr;
    filtmap=darklongdat.map(i).filtmap;
    wi=nfr*(nfr^2-1);
    [~,maski]=get_skymap(filtmap,ones(1024),4,5);
    [filtmap]=dc_offset_remove(filtmap,maski);
    dcstack=dcstack+(filtmap.*maski.*wi);
    maskstack=maskstack+maski.*wi;
end
dcstack=dcstack./maskstack;dcstack(find(dcstack~=dcstack))=0;
dcstack(find(dcstack==inf))=0;dcstack(find(dcstack==-inf))=0;
darklongdat.DCtemplate_nomask=dcstack;
%%
save(sprintf('%s/darklongdat',savedir),'darklongdat');
return