function get_chi2_fitg1diff(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using filtered(inst=1) or raw(inst=2)dark for RN, 
%fit ph with different G1.
% This code find the chi2 of flight vs RN+ph as a func of G1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixscale=7;
mypaths=get_paths(flight);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/noisemodel/diffG1fit/');
load(sprintf('%s/fwdat',loaddir),'fwdat');
load(sprintf('%s/noisemodel/diffCldat',loaddir),'diffCldat');

cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;

G1_arr=-4.5:0.1:-2.5;
nfr_arr=2:15;

chi2dat.G1_arr=G1_arr;
chi2dat.nfr_arr=nfr_arr;

for nfr=nfr_arr
    
%%% determine the field has that nfr
field_use=[];
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    if dt.nfrhalf>=nfr
        field_use=[field_use ifield];
    end
end
chi2dat.nfr(nfr).field_use=field_use;

%%%
chi2tot_arr=zeros(numel(field_use),numel(G1_arr));
subchi2tot_arr=zeros(numel(field_use),numel(G1_arr));
for ifieldcount=1:numel(field_use)
    ifield=field_use(ifieldcount);
    dt=get_dark_times(flight,inst,ifield);
    %%% get field info %%%
    loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
    load(strcat(loaddir,'bigmask'),'bigmask');
    load(strcat(loaddir,'flightmap'),'flightmap');

    %%% get fweight %%%
    fCl2d_arr=zeros(numel(dt.time),1024,1024);
    for i=1:numel(dt.time)
        load(strcat(loaddir,'labmap',num2str(i)),'labmap');
        filtmap=squeeze(labmap.filtmap_arr(nfr-1,:,:));
        [~,maskin1]=get_skymap(filtmap,bigmask,4,5);
        filtmap=filtmap-mean(filtmap(find(maskin1)));
        filtmap=filtmap.*maskin1;
        filtmap=dc_offset_remove(filtmap,maskin1).*maskin1;
        [~,~,~,~,~,~,fCl2d] = get_angular_spec(filtmap,filtmap,pixscale);
        fCl2d_arr(i,:,:)=fCl2d;
    end
    fCl2d_std=squeeze(std(fCl2d_arr));
    fw=(fftshift(fftshift(1./squeeze(fCl2d_std))))';
    
    %%% simulate normalized phCl(G1=-1) %%%
    slopemap=flightmap.rawmapf;
    [~,maskin1]=get_skymap(slopemap,bigmask,4,5);
    slopemean=mean(slopemap(find(maskin1)));
    phCl_arr=zeros(100,29);
    for i=1:100
        phmap1=photonnoise_realization(ones(1024).*slopemean,-1,nfr,frate);
        phmap2=photonnoise_realization(ones(1024).*slopemean,-1,nfr,frate);
        phmap=(phmap1-phmap2)./2;
        phmap=phmap-mean(phmap(find(maskin1)));phmap=phmap.*maskin1;
        [Cl,l] = get_angular_spec(phmap,phmap,pixscale,'w',fw);
        phCl_arr(i,:)=Cl;
    end
    phClavg=squeeze(mean(phCl_arr));
    phClstd=squeeze(std(phCl_arr));
    
    %%% get rnCl %%%
    rnCl_arr=squeeze(diffCldat.dark(ifield).wfCl_arr(:,nfr-1,:));
    rnClavg=squeeze(mean(rnCl_arr));
    rnClstd=squeeze(std(rnCl_arr));
    
    %%% get flight Cl %%%
    fCl=squeeze(diffCldat.flight(ifield).wfCl_arr(nfr-1,:));

    disp(sprintf('=============field%d,nfr=%d==============',ifield,nfr));
    %%% calculate chi2 %%%
    for iG1=1:numel(G1_arr)
        G1=G1_arr(iG1);
        nClavg=rnClavg+phClavg./(-G1);
        varnCl=rnClstd.^2+(phClstd./(-G1)).^2;
        chi2=(fCl-nClavg).^2./varnCl;

        chi2tot=sum(chi2(find(chi2==chi2)));
        subchi2tot=sum(chi2(end-4:end));
        disp(sprintf...
        ('G1=%.1f,chi2tot/dof=%.3e,subchi2tot/dof=%.3e',...
        G1,chi2tot/numel(find(find(chi2==chi2))),subchi2tot/5));
        chi2tot_arr(ifieldcount,iG1)=chi2tot;
        subchi2tot_arr(ifieldcount,iG1)=subchi2tot;
    end
end
chi2dat.nfr(nfr).chi2tot_arr=chi2tot_arr;
chi2dat.nfr(nfr).subchi2tot_arr=subchi2tot_arr;
end
save(strcat(savedir,'chi2dat'),'chi2dat');

return