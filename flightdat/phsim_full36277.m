%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This code run 100 sims of white noise photon noise map and get PS.
% - The mean of map is from alldat.dsubmean.
% - The calfactor G1 and G2 use the new calibrated values.
% - The ph sim map is clipped with the max and min of alldat.calmap5,
%   the clipping is not necessary at all. Only <10 pixels will be clip
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band=1;
beam=7;
iter_clip=3;
flight=36277;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;
G1=cp(band).apf2eps;G2=cp(band).eps2nWpm2ps;
frate=cp(band).framerate;%frame/sec

alldatdir='/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/';
load(strcat(alldatdir,'band',num2str(band),'_alldat'));
wdir=strcat('/Users/ytcheng/ciber/doc/20161102_36277/NoiseModel/TM',...
                num2str(band),'/full/PS2D/');
savedir=strcat('/Users/ytcheng/ciber/doc/20161102_36277/FlightDat/band',...
    num2str(band),'_ph/');
%%
for i=1:3
nframes=alldat(i).nfr;
Tint=nframes/frate; %sec

load(strcat(wdir,'b',num2str(band),'_i',...
      num2str(i), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';

calmap=alldat(i).calmap;
F=mean(calmap(find(calmap)))/G2;% mean of map (e-/s)
sig_photon=sqrt(6*F*(nframes^2+1)/5/Tint/(nframes^2-1));%e-/s

clipmax=(max(calmap(find(calmap~=0)))-mean(calmap(find(calmap))));
clipmin=(min(calmap(find(calmap~=0)))-mean(calmap(find(calmap))));

Mkk=alldat(i).wMkk;
mask=alldat(i).bigmask;

phCl_arr=zeros(100,29);
for j=1:100
disp(sprintf('field=%d,run=%d',i,j));

simmap=normrnd(0,sig_photon,1024);
simmap=simmap.*G2;
mask(simmap>clipmax | simmap<clipmin)=0;
cCl=get_Cl(simmap,mask,Mkk,beam,weight);
phCl_arr(j,:)=cCl;
end
save(strcat(savedir,'i',num2str(i),'phCl_arr'),'phCl_arr');
end
%% make plot
[~,l] = get_angular_spec(randn(1024),randn(1024),beam);

for i=1:3
figure
load(strcat(savedir,'i',num2str(i),'phCl_arr'),'phCl_arr');
for j=1:100
Cl=phCl_arr(j,:);
loglog(l,l.*(l+1).*Cl./2./pi,'k');hold on
end 
xlabel('$\ell$','interpreter','latex')
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW m^{-2} sr^{-1})$',...
        'interpreter','latex')
title(alldat(i).name);    
xlim([2e2,2e5]);ylim([1e-2,2e3]);
imname=strcat(savedir,'i',num2str(i),'ph_sim');
print(imname,'-dpng');
end