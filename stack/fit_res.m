%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% move the center or widen the PSF beam to fit the CB-UK residual map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
npix=400;
pixsize=0.7;
quad_arr=['A','B','C','D'];
mypaths=get_paths(flight);

psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/resmap/TM',...
    num2str(inst),'/');
savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/resmap/fitdat/TM',...
    num2str(inst),'/');

%%
load(strcat(savedir,'bestpar'),'bestpar');
for ifield=[4]%4:8
dt=get_dark_times(flight,inst,ifield);
for iquad=[3]%1:4
disp(sprintf('=================ifield=%d,iquad=%d=================',ifield,iquad))
quad=quad_arr(iquad);
resdat=fitsread(strcat(loaddir,dt.name,'_',quad,'_','resmap.fits'));

bestparam_norm=fitpsfdat(ifield).bestparam_norm;
A=bestparam_norm(1);
B=bestparam_norm(2);
sig=bestparam_norm(3);
r0=bestparam_norm(4);
alpha=bestparam_norm(5);
chi2best=bestparam_norm(6);

%%% 2MASS sim map PSF
radmap = make_radius_map(zeros(2*npix+1),npix,npix).*pixsize;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);

%%% CB PSF to be fit
fitrad=15;
sp=find(radmap<fitrad);
mask=zeros(size(resdat));
mask(sp)=1;
dof=numel(sp)-3;
bestpar(ifield).quad(iquad).chi2=1e10;

for rs=0.5:0.1:1.5
    disp(sprintf('rs=%.2f',rs));
    psfmapn = A*exp(-radmap.^2./2./(sig*rs)^2)+B./(1+(radmap./(r0*rs)).^alpha);
    norm1=sum(psfmapn(:));
    
    if ifield==5
        dxrange=4:0.1:7;%0:0.1:4;%1:0.1:2.7;
        dyrange=-4:0.1:0;%-4:0.1:-0.6;
    else
        dxrange=4:0.1:10;
        dyrange=-4:0.1:4;
    end
    
    for dx=dxrange
        for dy=dyrange
            radmap1 = make_radius_map(zeros(2*npix+1),npix-dy,npix-dx).*pixsize;
            psfmap1 = A*exp(-radmap1.^2./2./(sig*rs)^2)+...
                B./(1+(radmap1./(r0*rs)).^alpha);
            psfmap1=psfmap1./norm1;
            resmap=psfmap1-psfmap;
            resmod=resmap./psfmap(401,401);

            chi2=sum(sum((resdat-resmod).^2.*mask));
            
            if bestpar(ifield).quad(iquad).chi2>chi2
                bestpar(ifield).quad(iquad).dx=dx;
                bestpar(ifield).quad(iquad).dy=dy;
                bestpar(ifield).quad(iquad).rs=rs;
                bestpar(ifield).quad(iquad).chi2=chi2;
            end  
        end
    end
end
dx=bestpar(ifield).quad(iquad).dx;
dy=bestpar(ifield).quad(iquad).dy;
rs=bestpar(ifield).quad(iquad).rs;
chi2=bestpar(ifield).quad(iquad).chi2;
disp(sprintf('dx=%.1f,dy=%.1f,rs=%.2f,chi2/dof=%.2e',dx,dy,rs,chi2/dof))
end
end
save(strcat(savedir,'bestpar'),'bestpar');
%%
for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
for iquad=1:4
dx=bestpar(ifield).quad(iquad).dx;
dy=bestpar(ifield).quad(iquad).dy;
rs=bestpar(ifield).quad(iquad).rs;
chi2=bestpar(ifield).quad(iquad).chi2;
disp(sprintf('field%d,Q%d,dx=%.1f,dy=%.1f,rs=%.2f',ifield,iquad,dx,dy,rs))
   
end
disp('===========')
end
