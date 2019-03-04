flight=40030;
Nsims=100;% # of random position in the sims

npix=400;
quad_arr=['A','B','C','D'];


for inst=[1,2]
savedir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/forward_model/'));
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');
psfmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf/yt/inst',...
        num2str(inst),'/j0_14/');

for ifield=[8:-1:4]
dt = get_dark_times(flight,inst,ifield);

%%% get the CIBER stacked PSF map and profile %%%
radmap = make_radius_map(zeros(801),401,401).*0.7;
mapcb = zeros(801);
for quad=quad_arr
    q=fitsread(strcat(psfmapdir,dt.name,'_',quad,'.fits'));
    qh=fitsread(strcat(psfmapdir,dt.name,'_',quad,'hitmap.fits'));
    qmap=q./qh;
    bk=mean(qmap(find(radmap>100)));
    qmap=qmap-bk;
    mapcb=mapcb+qmap;
end
mapcb=mapcb./sum(mapcb(:));

profile = radial_prof(mapcb,ones(2*npix+1),401,401,1,15);
r_arr = profile.r.*0.7;
pcb_arr=(profile.prof)./profile.prof(1);
ecb_arr=profile.err./profile.prof(1);

%%% get the param for CIBER profile as prior %%%
bestparam = fitpsfdat(ifield).bestparam;
A_p=bestparam(1);
B_p=bestparam(2);
sig_p=bestparam(3);
r0_p=bestparam(4);
alpha_p=bestparam(5);

%%% iterate on A, B, alph %%%%
sig_arr = sig_p * (0.1:0.1:1);
r0_arr = r0_p * (0.1:0.1:1);
alpha_arr = alpha_p * (0.5:0.1:1.5);

sig_best = 0.;
r0_best = 0.;
alpha_best = 0.;
loss_best = 1e10;
stamp_best=zeros(2*npix+1,2*npix+1);
psim_best_arr = pcb_arr;
psfmap_best = mapcb;

% uniform between [npix+0.5,npix+0.5+10]
src_coord=npix+0.5+10*rand(2,Nsims);

for iscale=1:numel(sig_arr)
    sig = sig_arr(iscale);
    r0 = r0_arr(iscale);
for alpha=alpha_arr
    %%% get normalise amplitude A,B%%%
    radmap_cent = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*0.7;
    psfmap_cent = A_p*exp(-radmap_cent.^2./2./(sig)^2)+...
        B_p./(1+(radmap_cent./(r0)).^alpha);
    A=A_p/sum(psfmap_cent(:));
    B=B_p/sum(psfmap_cent(:));
    psfmap_cent = A*exp(-radmap_cent.^2./2./(sig)^2)+...
        B./(1+(radmap_cent./(r0)).^alpha);
    
    %%% do the random center sim for Nsims times %%%
    stamp=zeros(2*npix+1,2*npix+1);
    for isim = 1:Nsims
        xsrc=src_coord(1,isim);
        ysrc=src_coord(2,isim);
        radmap = make_radius_map(zeros(2*npix+10),xsrc,ysrc).*0.7;
        psfmap = A*exp(-radmap.^2./2./(sig)^2)+...
                 B./(1+(radmap./(r0)).^alpha);
        psfmap_coarse=rebin_map_coarse(psfmap,10);
        psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');

        stamp=stamp+psfmap_fine(round(xsrc)-npix:round(xsrc)+npix,...
                                round(ysrc)-npix:round(ysrc)+npix);
    end
    stamp=stamp./Nsims;
    
    %%% calculate the profile %%%
    profile = radial_prof(stamp,ones(2*npix+1),npix+1,npix+1,1,15);
    psim_arr=(profile.prof)./profile.prof(1);

    %%% calculate the chi2 %%%%
    sp = find(r_arr<5e1);
    loss = (log10(pcb_arr) - log10(psim_arr)).^2;
    loss = sum(loss(sp));
    if loss < loss_best
        loss_best = loss;
        sig_best = sig;
        r0_best = r0;
        alpha_best = alpha;
        stamp_best = stamp;
        psim_best_arr = psim_arr;
        psfmap_best = psfmap_cent;
    end
end
disp(sprintf('TM%d, field%d, sim %d/%d',inst,ifield,iscale,numel(sig_arr)))
end

disp(sprintf('TM%d, field%d',inst,ifield));

%%% save the best fit params
psfmodel(ifield).sig_best = sig_best;
psfmodel(ifield).r0_best = r0_best;
psfmodel(ifield).alpha_best = alpha_best;
psfmodel(ifield).r_arr = r_arr;
psfmodel(ifield).prof_best_arr = psim_best_arr;
psfmodel(ifield).psfmap = psfmap;

%%% plot the results
fig=figure;
setwinsize(gcf,1300,400)

v=[-1e-4,1e-4];
subplot(1,3,1)
x_arr=0.7*(-50:50);
imagesc(x_arr,x_arr,stamp_best(npix+1-50:npix+1+50,npix+1-50:npix+1+50)-...
    mapcb(npix+1-50:npix+1+50,npix+1-50:npix+1+50));
title('rebinned stack-CIBER stack');
xlabel('arcsec')
ylabel('arcsec')
caxis(v);
v=caxis;

subplot(1,3,2)
errorbar(r_arr,pcb_arr,ecb_arr,'ro-','Displayname','CIBER stack');hold on
semilogx(r_arr,psim_best_arr,'color',[0.00  0.50  0.00],'linewidth',2,...
    'Displayname','best fit model');

set(gca,'XScale','log','YScale','lin');
ylim([-0.1,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')
title(dt.name);
h=legend('show','Location','northeast');
set(h,'fontsize',10)
legend boxoff

subplot(1,3,3)
errorbar(r_arr,pcb_arr,ecb_arr,'ro-');hold on
errorbar(r_arr,-pcb_arr,ecb_arr,'bo-');hold on
plot(r_arr,psim_best_arr,'color',[0.00  0.50  0.00],'linewidth',2);
set(gca,'XScale','log','YScale','log');
ylim([1e-5,1.1])
xlim([4e-1,7e2])
xlabel('arcsec')
ylabel('<I_{stack}>')

drawnow
savename = strcat(savedir,'TM',num2str(inst),'_',dt.name,'_psfsim');
print(savename,'-dpng');


save(strcat(savedir,'TM',num2str(inst),'_psfmodel'),'psfmodel');

end

end

%% write the best fit param to txt file
flight=40030;
npix=400;
for inst=[1,2]
loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');
savedir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/forward_model/'));
load(strcat(savedir,'TM',num2str(inst),'_psfmodel'));
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    bestparam=fitpsfdat(ifield).bestparam;
    A_p = bestparam(1);
    B_p = bestparam(2);
    sig = psfmodel(ifield).sig_best;
    r0 = psfmodel(ifield).r0_best;
    alpha = psfmodel(ifield).alpha_best;
    
    % normalize A,B
    radmap_cent = make_radius_map(zeros(2*npix+1),npix+1,npix+1).*0.7;
    psfmap_cent = A_p*exp(-radmap_cent.^2./2./(sig)^2)+...
        B_p./(1+(radmap_cent./(r0)).^alpha);
    A=A_p/sum(psfmap_cent(:));
    B=B_p/sum(psfmap_cent(:));
    
    bestparam=[A B sig r0 alpha];
    dlmwrite(strcat(savedir,'TM',num2str(inst),'_',dt.name,'_bestparam.txt')...
    ,bestparam);
end
end
