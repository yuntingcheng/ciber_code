%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Get and save 12 frs 1st, 2nd, diff PS
% - Do the 1st, 2nd, diff for row and cal[(raw-DCtemplate)/FF]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
band=1;
beam=7;

halfdatdir='/Users/ytcheng/ciber/doc/20160920_FlightDiff/';
load(strcat(halfdatdir,'band',num2str(band),'_halfdat'),'halfdat');
weightdir=strcat('/Users/ytcheng/ciber/doc/20160905_NoiseModel/TM',...
            num2str(band),'/full/PS2D/');
rndir='/Users/ytcheng/ciber/doc/20160905_NoiseModel/diff/PS1D/';
savedir='/Users/ytcheng/ciber/doc/20160920_FlightDiff/';
%%

for i=1:8
%%% get weight
load(strcat(weightdir,'b',num2str(band),'_i',...
      num2str(i), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';
%%% get data
bigmask=halfdat(i).bigmask;
Mkk=halfdat(i).wMkk;
raw1=halfdat(i).raw1;raw2=halfdat(i).raw2;    
cal1=halfdat(i).cal1;cal2=halfdat(i).cal2;    
rawd=halfdat(i).rawd;cald=halfdat(i).cald;    

%%%
[r1Cl,~,r1Cl2d,l,binl]=get_Cl(raw1,bigmask,Mkk,beam,weight);
r1dCl=dCl_Knox(r1Cl,binl,beam);
[r2Cl,~,r2Cl2d]=get_Cl(raw2,bigmask,Mkk,beam,weight);
r2dCl=dCl_Knox(r2Cl,binl,beam);
[rdCl,~,rdCl2d]=get_Cl(rawd,bigmask,Mkk,beam,weight);
rddCl=dCl_Knox(rdCl,binl,beam);
[rdClnw]=get_Cl(rawd,bigmask,Mkk,beam,ones(1024));
rddClnw=dCl_Knox(rdClnw,binl,beam);

[c1Cl,~,c1Cl2d]=get_Cl(cal1,bigmask,Mkk,beam,weight);
c1dCl=dCl_Knox(c1Cl,binl,beam);
[c2Cl,~,c2Cl2d]=get_Cl(cal2,bigmask,Mkk,beam,weight);
c2dCl=dCl_Knox(c2Cl,binl,beam);
[cdCl,~,cdCl2d]=get_Cl(cald,bigmask,Mkk,beam,weight);
cddCl=dCl_Knox(cdCl,binl,beam);
[cdClnw]=get_Cl(cald,bigmask,Mkk,beam,ones(1024));
cddClnw=dCl_Knox(cdClnw,binl,beam);

%%% save all PS
psdat(i).l=l;psdat(i).name=halfdat(i).name;
psdat(i).r1Cl=r1Cl;psdat(i).r1dCl=r1dCl;psdat(i).r1Cl2d=r1Cl2d;
psdat(i).r2Cl=r2Cl;psdat(i).r2dCl=r2dCl;psdat(i).r2Cl2d=r2Cl2d;
psdat(i).rdCl=rdCl;psdat(i).rddCl=rddCl;psdat(i).rdCl2d=rdCl2d;
psdat(i).c1Cl=c1Cl;psdat(i).c1dCl=c1dCl;psdat(i).c1Cl2d=c1Cl2d;
psdat(i).c2Cl=c2Cl;psdat(i).c2dCl=c2dCl;psdat(i).c2Cl2d=c2Cl2d;
psdat(i).cdCl=cdCl;psdat(i).cddCl=cddCl;psdat(i).cdCl2d=cdCl2d;
psdat(i).rdClnw=rdClnw;psdat(i).rddClnw=rddClnw;
psdat(i).cdClnw=cdClnw;psdat(i).cddClnw=cddClnw;
end
save(strcat(savedir,'band',num2str(band),'_diffpsdat'),'psdat');
