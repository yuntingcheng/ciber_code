flight=40030;
inst=2;

if inst==1
    sin_freq=9.503;
elseif inst==2
    sin_freq=9.538;
end
mypaths=get_paths(flight);
darkdir=strcat(mypaths.Dark_Raw,'TM',num2str(inst),'/');
for ifield=1:8
dt=get_dark_times(flight,inst,ifield);
nfrhalf=dt.nfrhalf;
nfr_arr=2:dt.nfrhalf;

savedir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);
load(strcat(savedir,'maskin'),'maskin');

%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% filter the lab dark %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i=1;    
% read in raw frames   
disp(sprintf('retrive frames for ifield=%d,idark=%d/%d',...
    ifield,i,numel(dt.time)));
scanfile=dir(strcat(darkdir,'*',dt.time{i},'*'));
rawfr=zeros(dt.nfr,1024,1024);
for j=dt.frdown(i):dt.frup(i)
    fname=strcat(darkdir,scanfile(j).name);
    frame = imrotate(fitsread(fname),270);
    rawfr(j-dt.frdown(i)+1,:,:)=frame;
end

%%% line fit rawmap
[rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 
    
%%% filter lab dark
disp(sprintf('filter ifield=%d,idark=%d/%d',...
    ifield,i,numel(dt.time)));
[~,maskin1]=get_skymap(rawmap,maskin,2,5);

%%% do filtering
[filtfr] = imager_filtts(rawfr,'mapin',rawmap,'maskin',maskin1,...
 'offin',rawoff,'sin_freq',sin_freq,'verbose',1,'makeplot',1); 
  %}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% filter the flight %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[rawfr] = get_data_frames (inst,dt.name,'flight',flight,'verbose',0);
rawfr=rawfr(3:end,:,:);

disp(sprintf('linfit flight rawmap'));
[rawmap,rawoff]=linfit_map(rawfr,'verbose',0); 
[~,maskin1]=get_skymap(rawmap,maskin,2,5);

disp(sprintf('start filtering flight'));
[filtfr] = imager_filtts(rawfr,'mapin',rawmap,'maskin',maskin1,...
    'offin',rawoff,'sin_freq',sin_freq,'verbose',0,'makeplot',0);
sin_freq=9.518;
[filtfrnew] = imager_filtts(rawfr,'mapin',rawmap,'maskin',maskin1,...
    'offin',rawoff,'sin_freq',sin_freq,'verbose',0,'makeplot',0);

%%% line fit map
fr1=rawfr(1:nfrhalf,:,:);
fr2=rawfr(nfrhalf+1:2*nfrhalf,:,:);
[rawmap1] =linfit_map(fr1,'verbose',0);
[rawmap2] =linfit_map(fr2,'verbose',0);
[rawmapf] =linfit_map(rawfr,'verbose',0); 
flightmap.rawmap1=rawmap1;
flightmap.rawmap2=rawmap2;
flightmap.rawmapf=rawmapf;

fr1=filtfr(1:nfrhalf,:,:);
fr2=filtfr(nfrhalf+1:2*nfrhalf,:,:);
[filtmap1] =linfit_map(fr1,'verbose',0);
[filtmap2] =linfit_map(fr2,'verbose',0);
[filtmapf] =linfit_map(filtfr,'verbose',0);
flightmap.filtmap1=filtmap1;
flightmap.filtmap2=filtmap2;
flightmap.filtmapf=filtmapf;

fr1=filtfrnew(1:nfrhalf,:,:);
fr2=filtfrnew(nfrhalf+1:2*nfrhalf,:,:);
[filtmap1] =linfit_map(fr1,'verbose',0);
[filtmap2] =linfit_map(fr2,'verbose',0);
[filtmapf] =linfit_map(filtfr,'verbose',0);
flightmap.filtmap1new=filtmap1;
flightmap.filtmap2new=filtmap2;
flightmap.filtmapfnew=filtmapf;

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%s/maskdat',savedir),'maskdat');

rawdiff=(flightmap.rawmap1-flightmap.rawmap2)./2;
filtdiff=(flightmap.filtmap1-flightmap.filtmap2)./2;
filtdiffnew=(flightmap.filtmap1new-flightmap.filtmap2new)./2;
[~,maskin1]=get_skymap(filtdiff,maskdat.mask(ifield).bigmask,4,5);
[~,maskin1]=get_skymap(filtdiff,maskin1,4,5);

rawdiff=rawdiff-mean(rawdiff(find(maskin1)));rawdiff=rawdiff.*maskin1;
filtdiff=filtdiff-mean(filtdiff(find(maskin1)));filtdiff=filtdiff.*maskin1;
filtdiffnew=filtdiffnew-mean(filtdiffnew(find(maskin1)));
filtdiffnew=filtdiffnew.*maskin1;

[rCl,l,~,~,binl,~,rCl2d]=get_angular_spec(rawdiff,rawdiff,7);
[fCl,l,~,~,~,~,fCl2d]=get_angular_spec(filtdiff,filtdiff,7);
[fClnew,l,~,~,~,~,fCl2dnew]=get_angular_spec(filtdiffnew,filtdiffnew,7);
Cldat(ifield).rCl=rCl;
Cldat(ifield).fCl=rCl;
Cldat(ifield).fClnew=fClnew;
Cldat(ifield).rCl2d=rCl2d;
Cldat(ifield).fCl2d=rCl2d;
Cldat(ifield).fCl2dnew=fCl2dnew;
end
%%
figure
loglog(l,l.*(l+1).*rCl,'ko');hold on
loglog(l,l.*(l+1).*fCl,'b.');
loglog(l,l.*(l+1).*fClnew,'r+');
xlim([1e2,2e5]);
legend({'no filtered','filtered (dark freq)','filtered'},'location','southeast')
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell$',...
    'interpreter','latex','fontsize',18)
%%
ell=get_l(1024,1024,7,1);
fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
[x,y]=find(fmask);

figure
setwinsize(gcf,1000,300)
subplot(1,3,1)
imageclip(rCl2d);
v=caxis;
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('unfiltered')

subplot(1,3,2)
imageclip(fCl2d);
caxis(v);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('filter--9.538Hz(dark freq)')

subplot(1,3,3)
imageclip(fCl2dnew);
caxis(v);
xlim([min(x),max(x)]); ylim([min(y),max(y)]);
title('filter--9.518Hz')


%%
Cl2dtot=zeros(1024);
for ifield=4:8
ell=get_l(1024,1024,7,1);
fmask=zeros(1024);fmask((ell >= binl(22)) & (ell <= binl(22+1)))=1;
[x,y]=find(fmask);
Cl2dtot=Cl2dtot+Cldat(ifield).rCl2d;
figure
imageclip(Cldat(ifield).rCl2d(min(x):max(x),min(y):max(y)));hold on
colormap('parula')
[c,h]=contour(ell(min(x):max(x),min(y):max(y)),[binl(20),binl(21)],'color','m');
set(h,'linewidth',1.5);
[c,h]=contour(ell(min(x):max(x),min(y):max(y)),[binl(18),binl(19),binl(22),binl(23),binl(24)],'color','g');
set(h,'linewidth',1.5);

title(ifield);
savename=strcat('~/Desktop/field',num2str(ifield));
print(savename,'-dpng');close

end
%%

figure
imageclip(Cl2dtot(min(x):max(x),min(y):max(y)));hold on
colormap('parula')
[c,h]=contour(ell(min(x):max(x),min(y):max(y)),[binl(20),binl(21)],'color','m');
set(h,'linewidth',1.5);
[c,h]=contour(ell(min(x):max(x),min(y):max(y)),[binl(18),binl(19),binl(22),binl(23),binl(24)],'color','g');
set(h,'linewidth',1.5);

savename=strcat('~/Desktop/fieldstack');
print(savename,'-dpng');close
