%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit the spherical average PSF profile with an analytic function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
pixsize=0.7;
mypaths=get_paths(flight);

psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf/yt/inst',...
        num2str(inst),'/j0_14/');

savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
%% fit psf 1D profile
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
name=dt.name;

a=fitsread(strcat(psfdir,name,'_A.fits'));
b=fitsread(strcat(psfdir,name,'_B.fits'));
c=fitsread(strcat(psfdir,name,'_C.fits'));
d=fitsread(strcat(psfdir,name,'_D.fits'));
ah=fitsread(strcat(psfdir,name,'_Ahitmap.fits'));
bh=fitsread(strcat(psfdir,name,'_Bhitmap.fits'));
ch=fitsread(strcat(psfdir,name,'_Chitmap.fits'));
dh=fitsread(strcat(psfdir,name,'_Dhitmap.fits'));

psfdata=(a.*ah+b.*bh+c.*ch+d.*dh)./(ah+bh+ch+dh);

Starting(5)=max(psfdata(:));
[xindmax,yindmax]=find(psfdata==Starting(5));
Starting(3)=xindmax;Starting(4)=yindmax;
[spx,spy]=find(psfdata<1.1*exp(-0.5)*Starting(5)...
            & psfdata>0.9*exp(-0.5)*Starting(5));
spx=spx-xindmax;spy=spy-yindmax;
Starting(1)=mean(sqrt(spx(:).^2))
Starting(2)=mean(sqrt(spy(:).^2))

options=optimset('Display','iter');
gfitparam=fminsearch(@gaussfit2d,Starting,options,psfdata);
fitpsfdat(ifield).gfitparam=gfitparam;
end
save(strcat(savedir,'fitpsfdat'),'fitpsfdat');

%% write the best fit param to txt file for IDL read
load(strcat(savedir,'fitpsfdat'),'fitpsfdat');
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    
    gfitparam=fitpsfdat(ifield).gfitparam;    
    dlmwrite(strcat(savedir,dt.name,'_bestparam_g.txt'),gfitparam);
end