%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Get DC Tempalte by stacking couple sets of long integration time.
% - Get Inst Mask by masking out crazy pixels in DC Template.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
band=1;
%%
if flight==36277
darkdir=strcat('/Volumes/HD1TB/CIBER/data/Dark_Raw/362xx/TM',...
    num2str(band),'/');
savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/36277/';
elseif flight==40030
darkdir=strcat('/Volumes/HD1TB/CIBER/data/Dark_Raw/40030/TM',...
    num2str(band),'/');
savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/';
end
        
iter_clip=3;
sig_clip=3;
[time_arr,start_arr,end_arr]=get_dark_long(flight,band);

nfr_arr=end_arr-start_arr+1;
sig_arr=sqrt(1./nfr_arr./(nfr_arr.^2-1));
%% print out the quadrant mean
%{
for i=1:numel(time_arr)
time=time_arr{i};
scanfile=dir(strcat(darkdir,'*',time,'*'));

startfr=start_arr(i);endfr=end_arr(i);nfr=nfr_arr(i);
frames=zeros(nfr,1024,1024);
for j=startfr:endfr
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    frames(j-startfr+1,:,:)=frame;
end
darkmap=linfit_map(frames,'verbose',0);
%figure
%imageclip(darkmap);
%title(sprintf('%d',i));

mask=mask_clip(darkmap,iter_clip);

[m,m11,m12,m21,m22]=quadrant_mean(darkmap,mask);
pr=sprintf('%d,%s,%d-%d,m=%.2f,m11=%.2f,m12=%.2f,m21=%.2f,m22=%.2f',...
    i,time,startfr,endfr,m,m11,m12,m21,m22);
disp(pr);

end
%}
%% the index of bad dark, not using in DC template
if flight==36277
    switch band
    case 1
        badindex=[1,11,12,13,15,17];
    case 2
        badindex=[1,3,13,14,15,19,20];
    end

elseif flight==40030    
    switch band
    case 1
        badindex=[];
    case 2
        badindex=[9,11];
    end
end
%% stack DC map
use_arr=1:numel(time_arr);use_arr(badindex)=[];
dcstack=zeros(1024);
maskstack=zeros(1024);
for i=use_arr
time=time_arr{i};
scanfile=dir(strcat(darkdir,'*',time,'*'));

startfr=start_arr(i);endfr=end_arr(i);nfr=nfr_arr(i);
frames=zeros(nfr,1024,1024);
for j=startfr:endfr
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    frames(j-startfr+1,:,:)=frame;
end
darkmap=linfit_map(frames,'verbose',0);

mask=sigma_clip(darkmap,sig_clip,iter_clip);

figure
imageclip(darkmap.*mask);
title(sprintf('%d',i));

dcstack=dcstack+(darkmap.*mask./sig_arr(i)^2);
maskstack=maskstack+mask./sig_arr(i)^2;
pr=sprintf('%d,%s,%d-%d,nfr=%d',i,time,startfr,endfr,nfr);
disp(pr);

end

dcstack=dcstack./maskstack;dcstack(find(dcstack~=dcstack))=0;
dcstack(find(dcstack==inf))=0;dcstack(find(dcstack==-inf))=0;
%% make the mask_inst
mask=sigma_clip(dcstack,sig_clip,iter_clip);

sm_map=fillpadsmooth(dcstack,mask,2);
sm_mask=sigma_clip(sm_map,sig_clip,iter_clip);mask_inst=mask.*sm_mask;

%%% hand mask
mask_inst(513,:)=0;
mask_inst(:,1008:end)=0;
mask_inst(820:end,985:end)=0;
mask_inst(985:end,940:end)=0;
mask_inst(450:550,980:end)=0;

imageclip(dcstack.*mask_inst);
%%
DCtemplate=dcstack.*mask_inst;

imageclip(DCtemplate);
ylabel(colorbar,'ADU/fr','fontsize',18);
pltname=strcat(savedir,'band',num2str(band),'_DCtemplate');
print(pltname,'-dpng');%close

imageclip(mask_inst);
pltname=strcat(savedir,'band',num2str(band),'_mask_inst');
print(pltname,'-dpng');%close

save(strcat(savedir,'band',num2str(band),'_DCtemplate'),'DCtemplate');
save(strcat(savedir,'band',num2str(band),'_mask_inst'),'mask_inst');

