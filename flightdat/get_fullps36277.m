%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Get and save full frs auto PS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
band=2;
beam=7;
iter_clip=3;
flight=36277;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;

alldatdir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'));
wdir=strcat('/Users/ytcheng/ciber/doc/20161102_36277/NoiseModel/TM',...
                num2str(band),'/full/PS2D/');
rndir=strcat('/Users/ytcheng/ciber/doc/20161102_36277/NoiseModel/TM',...
                    num2str(band),'/full/PS1D/');
savedir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/';
%%
for i=1:3
bigmask=alldat(i).bigmask;
rawmap=alldat(i).rawmap.*bigmask.*cal;
calmap=alldat(i).calmap;
fieldname=alldat(i).name;
Mkk=alldat(i).wMkk;
%%%%!!!!!Mkkg makes negative Cl at BootesB, shoudl fiugre out....
Mkkg=alldat(i).wMkk;

%%% get weight
load(strcat(wdir,'b',num2str(band),'_i',...
      num2str(i), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';

%%% fit and sub plane in calmap
plane=plane_fit(calmap,bigmask);
scalmap=(calmap-plane).*bigmask;

%%% get 2D and weighted 1D PS
[rCl,~,rCl2d,l,binl]=get_Cl(rawmap,bigmask,Mkk,beam,weight);
drCl=dCl_Knox(rCl,binl,beam);
[cCl,~,cCl2d]=get_Cl(calmap,bigmask,Mkk,beam,weight);
dcCl=dCl_Knox(cCl,binl,beam);
[sCl,~,sCl2d]=get_Cl(scalmap,bigmask,Mkkg,beam,weight);
dsCl=dCl_Knox(sCl,binl,beam);
[sClnw]=get_Cl(scalmap,bigmask,Mkkg,beam,ones(1024));
dsClnw=dCl_Knox(sClnw,binl,beam);

%%% plot calmap and it's Cl2d
figure;
imageclip(scalmap);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_scalmap');
title(fieldname);print(imname,'-dpng');close

figure;
imageclip(sCl2d);
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_sCl2d');
title(fieldname);print(imname,'-dpng');close

figure;
imageclip(sCl2d.*weight');
imname=strcat(savedir,'band',num2str(band),'_plots/b',num2str(band),'_',...
                                    'i',num2str(i),'_sCl2d5w');
title(fieldname);print(imname,'-dpng');close

%%% save all PS
psdat(i).l=l;psdat(i).name=fieldname;
psdat(i).rCl=rCl;psdat(i).drCl=drCl;
psdat(i).cCl2d=cCl2d;
psdat(i).sCl2d=sCl2d;
psdat(i).cCl=cCl;psdat(i).dcCl=dcCl;%calmap PS / Knox err
psdat(i).sCl=sCl;psdat(i).dsCl=dsCl;%calmap-grad map PS/ knox err
psdat(i).sClnw=sClnw;psdat(i).dsClnw=dsClnw;%cal-grad (no Fourier W)
end
save(strcat(savedir,'band',num2str(band),'_fullpsdat'),'psdat');
