%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using filtered dark for RN, fit ph with different G1.
% This code find the chi2 of flight vs RN+ph as a func of G1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flight=40030;
inst=1;

cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;
G1_arr=-4.5:0.1:-2.5;
%%
slopemap_arr=zeros(8,1024,1024);
ifield_arr=[1,2,4,5,6,7,8];%%%%%%
for ifield=ifield_arr
    dt=get_dark_times(flight,inst,ifield);
    [flightfr] = get_data_frames...
        (inst,dt.name,'flight',flight,'verbose',0);
    flightfr=flightfr(3:end,:,:);
    slopemap=linfit_map(flightfr,'verbose',0);
    slopemap_arr(ifield,:,:)=slopemap;
end
%%
for nfr=18:19
%%% determine the field has that nfr
field_use=[];
for ifield=ifield_arr
    dt=get_dark_times(flight,inst,ifield);
    if dt.nfrhalf>=nfr
        field_use=[field_use ifield];
    end
end
chi2dat(nfr).field_use=field_use;

%%%
chi2tot_arr=zeros(numel(field_use),numel(G1_arr));
subchi2tot_arr=zeros(numel(field_use),numel(G1_arr));
for ifieldcount=1:numel(field_use)
    ifield=field_use(ifieldcount);
    dt=get_dark_times(flight,inst,ifield);
    %%% get field info %%%
    loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
           'sinfiltamp/DiffMap/field',num2str(ifield),'/');
    load(strcat(loaddir,'maskin'),'maskin');
    load(strcat(loaddir,'flightmap'),'flightmap');
    load(strcat(loaddir,'fCl2dstack_arr'),'fCl2dstack_arr');
    fCl2dstack=squeeze(fCl2dstack_arr(nfr-1,:,:));
    Wfilt=(fftshift(fftshift(1./fCl2dstack)))';

    %%% load rn and flight map %%%  
    rnmap_arr=zeros(numel(dt.time),1024,1024);
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labdat',num2str(i)),'labdat');
        rnmap=squeeze(labdat.filtmap_arr(nfr-1,:,:));
        rnmap_arr(i,:,:)=rnmap;
    end
    fmap=squeeze(flightmap.filtmap_arr(nfr-1,:,:));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    %%% get raw map for ph noise %%%
    dt=get_dark_times(flight,inst,ifield);
    [flightfr] = get_data_frames...
        (inst,dt.name,'flight',flight,'verbose',0);
    flightfr=flightfr(3:end,:,:);
    fr1=flightfr(1:nfr,:,:);
    fr2=flightfr(dt.nfrhalf+1:dt.nfrhalf+nfr,:,:);
    [slopemap1]=linfit_map(fr1,'verbose',0);
    [slopemap2]=linfit_map(fr2,'verbose',0);
    %}
    slopemap=squeeze(slopemap_arr(ifield,:,:));
    slopemean=mean(slopemap(find(maskin)));
    slopemap1=ones(1024)*slopemean;
    slopemap2=ones(1024)*slopemean;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% calculate chi2 %%%
    for iG1=1:numel(G1_arr)
        G1=G1_arr(iG1);
        [chi2]=chi2_g1fit_from_flgiht(rnmap_arr,fmap,slopemap1,...
            slopemap2,maskin,nfr,frate,G1,'w',Wfilt,'makeplot',0);
        chi2tot=sum(chi2(find(chi2==chi2)));
        subchi2tot=sum(chi2(end-4:end));
        disp(sprintf...
            ('field%d,G1=%.1f,nfr=%d,chi2tot=%.3e,subchi2tot=%.3e',...
            ifield,G1,nfr,chi2tot,subchi2tot));
        chi2tot_arr(ifieldcount,iG1)=chi2tot;
        subchi2tot_arr(ifieldcount,iG1)=subchi2tot;
    end
end
chi2dat(nfr).chi2tot_arr=chi2tot_arr;
chi2dat(nfr).subchi2tot_arr=subchi2tot_arr;
end
%%
savedir='/Users/ytcheng/ciber/doc/20170209_TsFilter/G1fit/';
%save(strcat(savedir,'chi2dat'),'chi2dat');
save(strcat(savedir,'chi2datw'),'chi2dat');