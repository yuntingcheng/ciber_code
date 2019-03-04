%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - This code run 100 sims of white noise ph noise map and get PS for diff.
% - The mean of map is (halfdat.mean1+halfdat.mean2)/2.
% - The calfactor G1 and G2 use the new calibrated values.
% - No clipping is applied as is has been verified not important at all
%   in phsim_full.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band=2;
beam=7;
iter_clip=3;
flight=40030;
cp=get_cal_params('flight',flight);
cal=cp(band).apf2eps.*cp(band).eps2nWpm2ps;
if band==1;G1=-2.6161;G2=cal./G1;end
if band==2;G1=-2.8296;G2=cal./G1;end
frate=cp(band).framerate;%frame/sec

weightdir='/Users/ytcheng/ciber/doc/20160905_NoiseModel/full/PS2D/';
halfdatdir='/Users/ytcheng/ciber/doc/20160920_FlightDiff/';
load(strcat(halfdatdir,'band',num2str(band),'_halfdat'),'halfdat');
savedir=strcat('/Users/ytcheng/ciber/doc/20160920_FlightDiff/band',...
    num2str(band),'_ph/');
%%
for i=1:8
nframes=halfdat(i).nfr;
Tint=nframes/frate; %sec

load(strcat(weightdir,'b',num2str(band),'_i',...
      num2str(i), '_std_noise'),'std_noise');
weight=(fftshift(fftshift(1./std_noise)))';

F=(halfdat(i).mean1+halfdat(i).mean2)./2./G2;% mean of map (e-/s)
sig_photon=sqrt(6*F*(nframes^2+1)/5/Tint/(nframes^2-1));%e-/s

Mkk=halfdat(i).wMkk;
mask=halfdat(i).bigmask;

phCl_arr=zeros(100,29);
for j=1:100
disp(sprintf('field=%d,run=%d',i,j));

simmap=normrnd(0,sig_photon,1024);
simmap=simmap.*G2;
cCl=get_Cl(simmap,mask,Mkk,beam,weight);
phCl_arr(j,:)=cCl;
end
save(strcat(savedir,'i',num2str(i),'phCl_arr'),'phCl_arr');
end
%% make plot
[~,l] = get_angular_spec(randn(1024),randn(1024),beam);

for i=1:8
figure
load(strcat(savedir,'i',num2str(i),'phCl_arr'),'phCl_arr');
for j=1:100
Cl=phCl_arr(j,:);
loglog(l,l.*(l+1).*Cl./2./pi,'k');hold on
end 
xlabel('$\ell$','interpreter','latex')
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW m^{-2} sr^{-1})$',...
        'interpreter','latex')
title(halfdat(i).name);    
xlim([2e2,2e5]);ylim([1e-1,1e4]);
imname=strcat(savedir,'i',num2str(i),'ph_sim');
print(imname,'-dpng');
end