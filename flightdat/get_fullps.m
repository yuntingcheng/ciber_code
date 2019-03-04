%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get and save full frs auto PS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
band=2;
beam=7;
iter_clip=3;
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'),'alldat');
weightdir='/Users/ytcheng/ciber/doc/20160905_NoiseModel/full/PS2D/';
rndir='/Users/ytcheng/ciber/doc/20160905_NoiseModel/full/PS1D/';
savedir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
%%
for i=1:8
bigmask=alldat(i).bigmask;
rawmap=alldat(i).rawmap.*bigmask.*cal;
calmap5=alldat(i).calmap5;
fieldname=alldat(i).name;
Mkk=alldat(i).wMkk;
%%%%!!!!!Mkkg makes negative Cl at BootesB, shoudl fiugre out....
Mkkg=alldat(i).wMkk;

%%% get weight
load(strcat(weightdir,'b',num2str(band),'_i',...
      num2str(i), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';

%%% fit and sub plane in calmap
plane5=plane_fit(calmap5,bigmask);
scalmap5=(calmap5-plane5).*bigmask;

%%% get 2D and weighted 1D PS
[rCl,~,rCl2d,l,binl]=get_Cl(rawmap,bigmask,Mkk,beam,weight);
drCl=dCl_Knox(rCl,binl,beam);
[cCl5,~,cCl2d5]=get_Cl(calmap5,bigmask,Mkk,beam,weight);
dcCl5=dCl_Knox(cCl5,binl,beam);
[sCl5,~,sCl2d5]=get_Cl(scalmap5,bigmask,Mkkg,beam,weight);
dsCl5=dCl_Knox(sCl5,binl,beam);
[sCl5nw]=get_Cl(scalmap5,bigmask,Mkkg,beam,ones(1024));
dsCl5nw=dCl_Knox(sCl5nw,binl,beam);

%%% plot calmap and it's Cl2d
figure;
imageclip(scalmap5);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_scalmap5');
title(fieldname);print(imname,'-dpng');close

figure;
imageclip(sCl2d5);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_sCl2d5');
title(fieldname);print(imname,'-dpng');close

figure;
imageclip(sCl2d5.*weight');
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_sCl2d5w');
title(fieldname);print(imname,'-dpng');close

%%% noise substraction
load(sprintf('%sb%d_i%d_wnClFclip_arr',rndir,band,i),'wnClFclip_arr');

nCl=mean(wnClFclip_arr);ndCl=std(wnClFclip_arr);
subcCl5=cCl5-nCl;subdcCl5=sqrt((dcCl5).^2+(ndCl).^2);
subsCl5=sCl5-nCl;subdsCl5=sqrt((dsCl5).^2+(ndCl).^2);

%%% save all PS
psdat(i).l=l;psdat(i).name=fieldname;
psdat(i).rCl=rCl;psdat(i).drCl=drCl;
% FF5 case
psdat(i).cCl2d5=cCl2d5;
psdat(i).sCl2d5=sCl2d5;
psdat(i).cCl5=cCl5;psdat(i).dcCl5=dcCl5;%calmap PS / Knox err
psdat(i).sCl5=sCl5;psdat(i).dsCl5=dsCl5;%calmap-grad map PS/ knox err
psdat(i).sCl5nw=sCl5nw;psdat(i).dsCl5nw=dsCl5nw;%cal-grad (no Fourier W)
psdat(i).subcCl5=subcCl5;psdat(i).subdcCl5=subdcCl5;% calmap-RN/err
psdat(i).subsCl5=subsCl5;psdat(i).subdsCl5=subdsCl5;% calmap-grad-RN/err
end
save(strcat(savedir,'band',num2str(band),'_fullpsdat'),'psdat');
