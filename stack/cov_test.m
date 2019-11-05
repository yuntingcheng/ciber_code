flight=40030;
inst=1;
mypaths=get_paths(flight);
pltsavedir = strcat(mypaths.alldat,'plots/cov_test/');
%% plot fit with fiducial value
ifield = 4;
masklim = false;
dx = 1200;
nbins = 25;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',loaddir,dt.name),'stackdatall');
load(sprintf('%s/excessdat_%s',loaddir,dt.name),'excessdatall');

figure
setwinsize(gcf,1200,500)
for im = 4%1:4
stackdat = stackdatall(im).stackdat;
excessdat = excessdatall(im).excessdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
rsub_arr = stackdat.rsub_arr;
r_arr = stackdat.r_arr;

% data excess
dat = stackdat.all.profcbgsub - stackdat.bg.profcbgsub;
edat = excessdat.excess.profcbgsub;
cov_mat = excessdat.excov.covcbsub;
cov_inv = inv(cov_mat);
covj_mat = excessdat.exjcov.covcbsub;
covj_inv = inv(covj_mat);
covji_mat = excessdat.exjicov.covcbsub;
covji_inv = inv(covji_mat);

% hsc clustering excess
weight = stackdat.all.profhitg;
hscclusdat = get_hsc_clus_prof(flight, inst, masklim, weight);
clus = hscclusdat(4).linefitmodelsub;

% PSF x window map
psfwin_arr = excessdat.psf.profcb;
radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfwin_map = spline(r_arr(1:18),psfwin_arr(1:18),radmap);
psfwin_map(radmap > r_arr(18)) = 0;
profile = radial_prof(psfwin_map,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwin_arr = profile.prof./profile.prof(1);
[psfwinsub_arr] = profile_radial_binning(psfwin_arr,weight,1);

% model map
rrmap = make_radius_map(zeros(201),101,101).*0.7;
[hsc_map,hsc_map_sub,R200,param_arr] = HSC_Wang19_prof(rrmap,im,true);

hsc_map1 = squeeze(hsc_map_sub(1,:,:));
psfwinhsc_map1 = conv2(psfwin_map,hsc_map1,'same');

hsc_map2 = squeeze(hsc_map_sub(2,:,:));
psfwinhsc_map2 = conv2(psfwin_map,hsc_map2,'same');

Ie1  = param_arr(1,1);
n1 = param_arr(1,2);
xe1 = param_arr(1,3);
Ie2  = param_arr(2,1);
n2 = param_arr(2,2);
xe2 = param_arr(2,3);
chi2best = 1e8;
subplot(2,4,im)
hsc_map2i = sersic(rrmap./R200,Ie2,n2,xe1+xe2);
psfwinhsc_map2i = conv2(psfwin_map, hsc_map2i,'same');
profile = radial_prof(psfwinhsc_map2i + psfwinhsc_map1,...
    ones(2*dx+1),dx+1,dx+1,1,nbins);
prof = profile.prof / profile.prof(1);
efull = prof.*dat(1) - psfwin_arr.*dat(1);
esub = profile_radial_binning(efull,weight,1);
clus_scale = 1;
m = clus_scale.*clus + esub;
diff = edat - m;
chi2 = diff*cov_inv*diff';
chi2j = diff*covj_inv*diff';
chi2ji = diff*covji_inv*diff';
chi2mat = (diff'*diff).*cov_inv;
fprintf('xe = %.3f, clus_scale = %d,chi2 = %.3e,chi2j = %.3e,chi2ji = %.3e\n',...
    xe2,clus_scale,chi2,chi2j,chi2ji); 
loglog(rsub_arr,m);hold on
text(1e2,3e1,sprintf('chi^2_0 = %.1f \n chi^2_1 = %.1f \n chi^2_2 = %.1f',...
    chi2,chi2j,chi2ji));
text(4e-1,1e-1,sprintf('xe2 = %.4f \n A_{clus} = %d',xe2,clus_scale));

err = sqrt(diag(cov_mat))';
errl = err;
errl(errl > abs(edat)) = abs(edat(errl > abs(edat))) - 1e-10;
loglog(rsub_arr,edat,'k.','markersize',10);hold on
loglog(rsub_arr,-edat,'ko','markersize',5);
errorbar(rsub_arr,edat,errl,err,'k.','markersize',10);
errorbar(rsub_arr,-edat,errl,err,'ko','markersize',5);
xlim([4e-1,1.1e3])
ylim([1e-2,3e2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('Excess I [nW/m^2/sr]', 'fontsize',12);
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

subplot(2,4,im+4)
imageclip(chi2mat);
xticks([3:4:15])
xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
xtickangle(45)
yticks([3:4:15])
yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
ytickangle(45)
title(strcat(' \Delta_i C^{-1}_{ij} \Delta_j'));
end
savename = strcat(pltsavedir,'fit_fiducial'); 
print(savename,'-dpng');close
%% fit param
ifield = 4;
masklim = false;
dx = 1200;
nbins = 25;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',loaddir,dt.name),'stackdatall');
load(sprintf('%s/excessdat_%s',loaddir,dt.name),'excessdatall');
figure
setwinsize(gcf,1200,500)

for im = 1:4
stackdat = stackdatall(im).stackdat;
excessdat = excessdatall(im).excessdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
rsub_arr = stackdat.rsub_arr;
r_arr = stackdat.r_arr;

% data excess
dat = stackdat.all.profcbgsub - stackdat.bg.profcbgsub;
edat = excessdat.excess.profcbgsub;
cov_mat = excessdat.excov.covcbsub;
cov_inv = inv(cov_mat);
covj_mat = excessdat.exjcov.covcbsub;
covj_inv = inv(covj_mat);
covji_mat = excessdat.exjicov.covcbsub;
covji_inv = inv(covji_mat);

% hsc clustering excess
weight = stackdatall(im).stackdat.all.profhitg;
hscclusdat = get_hsc_clus_prof(flight, inst, masklim, weight);
clus = hscclusdat(4).linefitmodelsub;
clusfull = hscclusdat(4).linefitmodel;

% PSF x window map
psfwin_arr = excessdat.psf.profcb;
radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfwin_map = spline(r_arr(1:18),psfwin_arr(1:18),radmap);
psfwin_map(radmap > r_arr(18)) = 0;
profile = radial_prof(psfwin_map,ones(2*dx+1),dx+1,dx+1,1,nbins);
psfwin_arr = profile.prof./profile.prof(1);
[psfwinsub_arr] = profile_radial_binning(psfwin_arr,weight,1);

% model map
rrmap = make_radius_map(zeros(201),101,101).*0.7;
[hsc_map,hsc_map_sub,R200,param_arr] = HSC_Wang19_prof(rrmap,im,true);
% clus(rsub_arr<R200)=0;
% clusfull(r_arr<R200)=0;
clus(rsub_arr<10)=0;
clusfull(r_arr<10)=0;

hsc_map1 = squeeze(hsc_map_sub(1,:,:));
psfwinhsc_map1 = conv2(psfwin_map,hsc_map1,'same');

hsc_map2 = squeeze(hsc_map_sub(2,:,:));
psfwinhsc_map2 = conv2(psfwin_map,hsc_map2,'same');

Ie1  = param_arr(1,1);
n1 = param_arr(1,2);
xe1 = param_arr(1,3);
Ie2  = param_arr(2,1);
n2 = param_arr(2,2);
xe2 = param_arr(2,3);
xe2_arr = [0.005,0.01,0.04,0.07,0.10];
clus_scale_arr = [1,3,5,7];
chi2_arr = zeros([numel(xe2_arr,clus_scale_arr)]);
subplot(2,4,im)
bestfitdat(im).chi2 = inf;
bestfitdat(im).cov = cov_mat;
bestfitdat(im).covj = covj_mat;
bestfitdat(im).covji = covji_mat;

for ixe2= 1:numel(xe2_arr)
    xe2 = xe2_arr(ixe2);
    hsc_map2i = sersic(rrmap./R200,Ie2,n2,xe1+xe2);
    psfwinhsc_map2i = conv2(psfwin_map, hsc_map2i,'same');
    profile = radial_prof(psfwinhsc_map2i + psfwinhsc_map1,...
        ones(2*dx+1),dx+1,dx+1,1,nbins);
    prof = profile.prof / profile.prof(1);
    efull = prof.*dat(1) - psfwin_arr.*dat(1);
    esub = profile_radial_binning(efull,weight,1);
    for iclus_scale = 1:numel(clus_scale_arr)
        clus_scale = clus_scale_arr(iclus_scale);
        m = clus_scale.*clus + esub;
        diff = edat - m;
        chi2 = diff*cov_inv*diff';
        chi2_arr(ixe2,iclus_scale) = chi2;
        fprintf('xe = %.4f, clus_scale = %d,chi2 = %.2e\n',...
            xe2,clus_scale,chi2);
        if xe2==0.005
            loglog(rsub_arr,clus_scale.*clus,'DisplayName',...
            sprintf('A_{clus} = %d',clus_scale));hold on
        end
        
    if chi2<bestfitdat(im).chi2
        bestfitdat(im).chi2 = chi2;
        bestfitdat(im).xe2 = xe2;
        bestfitdat(im).clus_scale = clus_scale;
        bestfitdat(im).gal_model = esub;
        bestfitdat(im).clus_model = clus_scale.*clus;
        bestfitdat(im).model = clus_scale.*clus + esub;
        bestfitdat(im).gal_model_full = efull;
        bestfitdat(im).clus_model_full = clus_scale.*clusfull;
        bestfitdat(im).model_full = clus_scale.*clusfull + efull;
        bestfitdat(im).dat = edat;
        bestfitdat(im).daterr = sqrt(diag(cov_mat))';
    end
        
    end
    
    loglog(rsub_arr,esub,'DisplayName',...
    sprintf('xe2 = %.2f',xe2));hold on
end

p1=loglog(rsub_arr, bestfitdat(im).model, 'k', ...
    'linewidth',3 ,'DisplayName','best fit');
p1.Color(4) = 0.3;

if im==4
    h=legend('show','Location','northeast');
    set(h,'fontsize',10)
    legend boxoff
end

err = sqrt(diag(cov_mat))';
errl = err;
errl(errl > abs(edat)) = abs(edat(errl > abs(edat))) - 1e-10;
loglog(rsub_arr,edat,'k.','markersize',10);hold on
loglog(rsub_arr,-edat,'ko','markersize',5);
errorbar(rsub_arr,edat,errl,err,'k.','markersize',10);
errorbar(rsub_arr,-edat,errl,err,'ko','markersize',5);
xlim([4e-1,1.1e3])
ylim([1e-2,1e2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('Excess I [nW/m^2/sr]', 'fontsize',12);
title(strcat(num2str(m_min),'<m<',num2str(m_max)),'fontsize',15);

subplot(2,4,im+4)
for iclus_scale=1:numel(clus_scale_arr)
    clus_scale = clus_scale_arr(iclus_scale);
    plot(xe2_arr, chi2_arr(:,iclus_scale),'o-','DisplayName',...
            sprintf('A_{clus} = %d',clus_scale));hold on
end



xlabel('xe2', 'fontsize',15);
ylabel('\chi^2', 'fontsize',12);
title(strcat('\chi^2'));
if im==4
    h=legend('show','Location','northwest');
    set(h,'fontsize',10)
    legend boxoff
end

end
savename = strcat(pltsavedir,'fit_params'); 
print(savename,'-dpng');close
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check cov from 3 jackknife scheme are consistent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for im=3%1:4
    figure
    setwinsize(gcf,1200,500)
    
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    rsub_arr = stackdat.rsub_arr;
    r_arr = stackdat.r_arr;

    cov_mat = excessdat.excov.covcbsub;
    cov_inv = inv(cov_mat);
    covj_mat = excessdat.exjcov.covcbsub;
    covj_inv = inv(covj_mat);
    covji_mat = excessdat.exjicov.covcbsub;
    covji_inv = inv(covji_mat);
    
    subplot(2,4,1)
    edat = excessdat.excess.profcbg;
    edat_err = sqrt(diag(excessdat.excov.covcb))';
    m = bestfitdat(im).model_full;
    loglog(r_arr,m,'k');hold on
    loglog(r_arr,edat,'k.','markersize',10);
    loglog(r_arr,-edat,'ko','markersize',5);
    edat_errl = edat_err;
    edat_errl(edat_errl>abs(edat)) = abs(edat(edat_errl > abs(edat)))-1e-10;
    errorbar(r_arr,edat,edat_errl,edat_err,'k.','markersize',10);
    errorbar(r_arr,-edat,edat_errl,edat_err,'ko','markersize',5);    
    xlim([4e-1,1e3])
    ylim([1e-2,1e2])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('Excess I [nW/m^2/sr]', 'fontsize',12);
    title('data, best fit (25 bins)','fontsize',15);

    subplot(2,4,5)
    edat = excessdat.excess.profcbgsub;
    edat_err = sqrt(diag(excessdat.excov.covcbsub))';
    edat = bestfitdat(im).dat;
    edat_err = bestfitdat(im).daterr;
    m = bestfitdat(im).model;
    diff = edat - m;
    loglog(rsub_arr,m,'k');hold on
    loglog(rsub_arr,edat,'k.','markersize',10);
    loglog(rsub_arr,-edat,'ko','markersize',5);
    edat_errl = edat_err;
    edat_errl(edat_errl>abs(edat)) = abs(edat(edat_errl > abs(edat)))-1e-10;
    errorbar(rsub_arr,edat,edat_errl,edat_err,'k.','markersize',10);
    errorbar(rsub_arr,-edat,edat_errl,edat_err,'ko','markersize',5);    
    xlim([4e-1,1.1e3])
    ylim([1e-2,1e2])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('Excess I [nW/m^2/sr]', 'fontsize',12);
    title('data, best fit (15 bins)','fontsize',15);
       
    for itype=1:3  
        switch itype
        case 1
            exdat = excessdat.datcov;
            name = 'data';
        case 2
            exdat = excessdat.bgcov;
            name = 'background';
        case 3
            exdat = excessdat.psfcov;
            name = 'PSF';
        end

        subplot(2,4,itype+5)
        imageclip(exdat.covcbsub);
        title(strcat('Cov--',name,' (15 bins)'),'fontsize',15);
        xticks([3:4:15])
        xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        xtickangle(45)
        yticks([3:4:15])
        yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        ytickangle(45)        
        if itype==1
            v1 = caxis;
        else
            caxis(v1);
        end
        
        subplot(2,4,itype+1)
        imageclip(exdat.covcb);
        title(strcat('Cov--',name,' (25 bins)'),'fontsize',15);
        xticks([6:6:25])
        xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
            num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
        xtickangle(45)
        yticks([6:6:25])
        yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
            num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
        ytickangle(45)  
        caxis(v1);        
        
    end
    suptitle(strcat(num2str(m_min),'<m<',num2str(m_max)));
    savename = sprintf('%s%s_%d_%d_cov',pltsavedir,dt.name,m_min,m_max);
    print(savename,'-dpng');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    for itype=1:3  
        switch itype
        case 1
            exdat = excessdat.excov;
            name = 'exess';
        case 2
            exdat = excessdat.exjcov;
            name = 'excess jackknife';
        case 3
            exdat = excessdat.exjicov;
            name = 'excess jackknife w/ norm.';
        end
        
        figure
        setwinsize(gcf,900,500)
        subplot(2,3,2)
        imageclip(exdat.covcbsub);
        title(strcat('Cov (15 bins)'),'fontsize',15);
        xticks([3:4:15])
        xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        xtickangle(45)
        yticks([3:4:15])
        yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        ytickangle(45)        
        if itype==1
            v1 = caxis;
        else
            caxis(v1);
        end
        
        subplot(2,3,1)
        imageclip(exdat.covcb);
        title(strcat('Cov (25 bins)'),'fontsize',15);
        xticks([6:6:25])
        xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
            num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
        xtickangle(45)
        yticks([6:6:25])
        yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
            num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
        ytickangle(45)  
        caxis(v1);
        
        subplot(2,3,5)
        imageclip(normalize_cov(exdat.covcbsub));
        title(strcat('Norm. Cov (15 bins)'),'fontsize',15);
        xticks([3:4:15])
        xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        xtickangle(45)
        yticks([3:4:15])
        yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        ytickangle(45)        
        caxis([-1,1]);
        
        subplot(2,3,4)
        imageclip(normalize_cov(exdat.covcb));
        title(strcat('Norm. Cov (25 bins)'),'fontsize',15);
        xticks([6:6:25])
        xticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
            num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
        xtickangle(45)
        yticks([6:6:25])
        yticklabels({num2str(r_arr(6),'%.1e'),num2str(r_arr(12),'%.1e'),...
            num2str(r_arr(18),'%.1e'),num2str(r_arr(24),'%.1e')});
        ytickangle(45)  
        caxis([-1,1]);

        subplot(2,3,3)
        imageclip(inv(exdat.covcbsub));
        title(strcat('C^{-1}_{ij}'),'fontsize',15);        
        xticks([3:4:15])
        xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        xtickangle(45)
        yticks([3:4:15])
        yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        ytickangle(45)
        if itype==1
            v2 = caxis;
        else
            caxis(v2);
        end
        
        subplot(2,3,6)
        imageclip((diff'*diff).*inv(exdat.covcbsub));
        title(strcat('\Delta_i C^{-1}_{ij} \Delta_j'),'fontsize',15);                
        xticks([3:4:15])
        xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        xtickangle(45)
        yticks([3:4:15])
        yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
            num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
        ytickangle(45)
        if itype==1
            v3 = caxis;
        else
            caxis(v3);
        end
       suptitle(strcat(num2str(m_min),'<m<',num2str(m_max),', ',name));
       savename = sprintf('%s%s_%d_%d_excov%d',...
           pltsavedir,dt.name,m_min,m_max,itype-1);
       print(savename,'-dpng');

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%    
    chi2mat = (diff'*diff).*inv(excessdat.excov.covcbsub);
    chi2matj = (diff'*diff).*inv(excessdat.exjcov.covcbsub);
    chi2matji = (diff'*diff).*inv(excessdat.exjicov.covcbsub);
    figure
    semilogx(rsub_arr,sum(chi2mat),'k');hold on
    semilogx(rsub_arr,sum(chi2matj),'b');
    semilogx(rsub_arr,sum(chi2matji),'r');
    h=legend({'\chi^2','\chi^2(jackknife)','\chi^2(jackknife w/ norm)'},...
        'Location','northeast');
    hline(sum(chi2mat(:))/13,'k--');
    hline(sum(chi2mat(:))/13,'b--');
    hline(sum(chi2mat(:))/13,'r--');
    set(h,'fontsize',8)
    xlim([4e-1,1.1e3])
    ylim([-1,7])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('\Sigma_j( \Delta_i C^{-1}_{ij} \Delta_j)', 'fontsize',12);
    savename = sprintf('%s%s_%d_%d_chi2cols',pltsavedir,dt.name,m_min,m_max);
    print(savename,'-dpng');

end
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cov are stable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for im=3%1:4
    figure
    setwinsize(gcf,1200,500)
    
    stackdat = stackdatall(im).stackdat;
    excessdat = excessdatall(im).excessdat;
    m_min = stackdat.m_min;
    m_max = stackdat.m_max;
    rsub_arr = stackdat.rsub_arr;
    r_arr = stackdat.r_arr;

    C = excessdat.excov.covcbsub;
    Cnorm = normalize_cov(C);
    Cinv = inv(C);
    d = excessdat.excess.profcbgsub;
    d_err = sqrt(diag(excessdat.excov.covcbsub))';
    m = bestfitdat(im).model;
    diff = d - m;
    
    %%% test invert subregion
    Cinv1 = inv(C(1:5,1:5));
    Cinv2 = inv(C(6:10,6:10));
    Cinv3 = inv(C(11:15,11:15));
    Cinvs = zeros(15);
    Cinvs(1:5,1:5) = inv(C(1:5,1:5));
    Cinvs(6:10,6:10) = inv(C(6:10,6:10));
    Cinvs(11:15,11:15) = inv(C(11:15,11:15));
    
    chi2mat = (diff'*diff).*Cinv;
    chi2mats = (diff'*diff).*Cinvs;
    
    setwinsize(gcf,900,500)
    subplot(2,3,1)
    imageclip(Cnorm);
    title('Norm. Cov','fontsize',15);
    caxis([-1,1]);
    vline([5.5,10.5],'r--');
    hline([5.5,10.5],'r--');
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    
    subplot(2,3,2)
    imageclip(Cinv);
    title(strcat('C^{-1}_{ij}'),'fontsize',15);
    v=caxis;
    vline([5.5,10.5],'r--');
    hline([5.5,10.5],'r--');
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,3,3)
    imageclip(Cinvs);
    title(strcat('C^{-1}_{ij} (subregion)'),'fontsize',15);
    caxis(v);
    vline([5.5,10.5],'r--');
    hline([5.5,10.5],'r--');
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    
    subplot(2,3,5)
    imageclip(chi2mat);
    title(strcat('\Delta_i C^{-1}_{ij} \Delta_j'),'fontsize',15);                
    v=caxis;
    vline([5.5,10.5],'r--');
    hline([5.5,10.5],'r--');
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,3,6)
    imageclip(chi2mats);
    title(strcat('\Delta_i C^{-1}_{ij} \Delta_j (subregion)'),'fontsize',15);        
    caxis(v);
    vline([5.5,10.5],'r--');
    hline([5.5,10.5],'r--');
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    
    subplot(2,3,4)
    semilogx(rsub_arr,sum(chi2mat),'k');hold on
    semilogx(rsub_arr,sum(chi2mats),'r');
    h=legend({'\chi^2','\chi^2(subregion)'},'Location','northeast');
    hline(sum(chi2mat(:))/13,'k--');
    hline(sum(chi2mats(:))/13,'r--');
    set(h,'fontsize',8)
    xlim([4e-1,1.1e3])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('\Sigma_j( \Delta_i C^{-1}_{ij} \Delta_j)', 'fontsize',12);
    
    savename = strcat(pltsavedir,'cov_subreg'); 
    print(savename,'-dpng');close

    %%%% test adding noise to C and invert %%%
    Cn1 = C;
    Cn1(12,1:4) = 0;
    Cn1(1:4,12) = 0;
    chi2matn1 = (diff'*diff).*inv(Cn1);
    Cn2 = C;
    Cn2(1:4,6:7) = -abs(Cn2(1:4,6:7));
    Cn2(6:7,1:4) = -abs(Cn2(6:7,1:4));
    chi2matn2 = (diff'*diff).*inv(Cn2);
    
    figure
    setwinsize(gcf,1200,500)
    subplot(2,4,1)
    imageclip(Cnorm);
    title('Norm. Cov','fontsize',15);
    caxis([-1,1]);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,4,5)
    imageclip(Cinv);
    title('C^{-1}','fontsize',15);
    v=caxis;
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    
    subplot(2,4,2);
    imageclip(normalize_cov(Cn1));hold on
    title('Norm. Cov 1','fontsize',15);
    plot([0.5,4.5],[11.5,11.5],'r');
    plot([0.5,4.5],[12.5,12.5],'r');
    plot([4.5,4.5],[11.5,12.5],'r');
    plot([11.5,11.5],[0.5,4.5],'r');
    plot([12.5,12.5],[0.5,4.5],'r');
    plot([11.5,12.5],[4.5,4.5],'r');
    caxis([-1,1]);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,4,6);
    imageclip(inv(Cn1));hold on
    plot([0.5,4.5],[11.5,11.5],'r');
    plot([0.5,4.5],[12.5,12.5],'r');
    plot([4.5,4.5],[11.5,12.5],'r');
    plot([11.5,11.5],[0.5,4.5],'r');
    plot([12.5,12.5],[0.5,4.5],'r');
    plot([11.5,12.5],[4.5,4.5],'r');
    caxis(v);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
   
    subplot(2,4,3);
    imageclip(normalize_cov(Cn2));hold on
    title('Norm. Cov 2','fontsize',15);
    plot([0.5,4.5],[5.5,5.5],'r');
    plot([0.5,4.5],[7.5,7.5],'r');
    plot([4.5,4.5],[5.5,7.5],'r');
    plot([5.5,5.5],[0.5,4.5],'r');
    plot([7.5,7.5],[0.5,4.5],'r');
    plot([5.5,7.5],[4.5,4.5],'r');
    caxis([-1,1]);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,4,7);
    imageclip(inv(Cn2));hold on
    plot([0.5,4.5],[5.5,5.5],'r');
    plot([0.5,4.5],[7.5,7.5],'r');
    plot([4.5,4.5],[5.5,7.5],'r');
    plot([5.5,5.5],[0.5,4.5],'r');
    plot([7.5,7.5],[0.5,4.5],'r');
    plot([5.5,7.5],[4.5,4.5],'r');
    caxis(v);
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)
    
    subplot(2,4,8)
    semilogx(rsub_arr,sum(chi2mat),'k');hold on
    semilogx(rsub_arr,sum(chi2matn1),'b');
    semilogx(rsub_arr,sum(chi2matn2),'r');
    h=legend({'\chi^2','\chi^2_1','\chi^2_2'},'Location','northeast');
    hline(sum(chi2mat(:))/13,'k--');
    hline(sum(chi2matn1(:))/13,'b--');
    hline(sum(chi2matn2(:))/13,'r--');
    set(h,'fontsize',8)
    xlim([4e-1,1.1e3])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('\Sigma_j( \Delta_i C^{-1}_{ij} \Delta_j)', 'fontsize',12);

    savename = strcat(pltsavedir,'cov_noise'); 
    print(savename,'-dpng');close    
    %%%%% rebinning test %%%%
    
    chi2mat = (diff'*diff).*Cinv;
    w = zeros([1,15]);
    w(2:14) = stackdat.all.profhitg(7:19);
    w(1) = sum(stackdat.all.profhitg(1:5));
    w(15) = sum(stackdat.all.profhitg(20:25));
    
    [diffs1,Cs1] = profile_cov_binning(diff,C,w,1:5);
    [rsub_arr1,~] = profile_cov_binning(rsub_arr,C,w,1:5);
    chi2mat1 = (diffs1'*diffs1).*inv(Cs1);

    [diffs2,Cs2] = profile_cov_binning(diff,C,w,11:15);
    [rsub_arr2,~] = profile_cov_binning(rsub_arr,C,w,11:15);
    chi2mat2 = (diffs2'*diffs2).*inv(Cs2);

    figure
    setwinsize(gcf,1200,500)
    subplot(2,4,1)
    imageclip(normalize_cov(C));
    title('Norm. Cov','fontsize',15);
    caxis([-1,1])
    subplot(2,4,2)
    imageclip(normalize_cov(Cs1));
    title('Norm. Cov (bin 1-5)','fontsize',15);
    caxis([-1,1])
    subplot(2,4,3)
    imageclip(normalize_cov(Cs2));
    title('Norm. Cov (bin 11-15)','fontsize',15);
    caxis([-1,1])
    subplot(2,4,5)
    imageclip(inv(C));
    title('Inv Cov','fontsize',15);
    v=caxis;
    subplot(2,4,6)
    imageclip(inv(Cs1));
    title('Inv Cov (bin 1-5)','fontsize',15);
    caxis(v)
    subplot(2,4,7)
    imageclip(inv(Cs2));
    title('Inv Cov (bin 11-15)','fontsize',15);
    caxis(v)
    
    
    subplot(2,4,8)
    semilogx(rsub_arr,sum(chi2mat),'k.-');hold on
    semilogx(rsub_arr1,sum(chi2mat1),'b.-');
    semilogx(rsub_arr2,sum(chi2mat2),'r.-');
    h=legend({'\chi^2','\chi^2 (bin 1-5)','\chi^2_2 (bin 11-15)'},...
        'Location','northwest');
    hline(sum(chi2mat(:))/13,'k--');
    hline(sum(chi2mat1(:))/13,'b--');
    hline(sum(chi2mat2(:))/13,'r--');
    set(h,'fontsize',8)
    xlim([4e-1,1.1e3])
    xlabel('r [arcsec]', 'fontsize',15);
    ylabel('\Sigma_j( \Delta_i C^{-1}_{ij} \Delta_j)', 'fontsize',12);
    
    savename = strcat(pltsavedir,'cov_rebin'); 
    print(savename,'-dpng');close

end
%% Two Toy Model Cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% toy model chi2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flight=40030;
inst=1;
mypaths=get_paths(flight);
pltsavedir = strcat(mypaths.alldat,'plots/cov_test/');
ifield = 4;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',loaddir,dt.name),'stackdatall');
load(sprintf('%s/excessdat_%s',loaddir,dt.name),'excessdatall');

figure
setwinsize(gcf,1400,300)
im = 2;

stackdat = stackdatall(im).stackdat;
excessdat = excessdatall(im).excessdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
rsub_arr = stackdat.rsub_arr;

cov = excessdat.excov.covcbsub;
covi = inv(cov);
d = excessdat.excess.profcbgsub;
e = sqrt(diag(cov))';

subplot(1,4,1)
loglog(rsub_arr, d, 'k.', 'markersize', 10);hold on
loglog(rsub_arr, e,'b.-');
loglog(rsub_arr, sqrt(1./diag(covi)),'r.-');
h=legend({'Excess','sqrt(diag(Excess Cov))','1/sqrt(diag(Excess Cov^{-1}))'},...
        'Location','northeast');
set(h,'fontsize',10)
loglog(rsub_arr, -d, 'ko', 'markersize', 5);hold on
xlim([4e-1,1.1e3])
ylim([1e-1,1e2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('I [nW/m^2/sr]', 'fontsize',15);


diff_list = [];

diff = e;
diff_list = [diff_list;diff];

diff = e;
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

diff = e.*(0.1:0.1:1.5);
diff_list = [diff_list;diff];

diff = e.*(0.1:0.1:1.5);
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

diff = e.*(1.5:-0.1:0.1);
diff_list = [diff_list;diff];

diff = e.*(1.5:-0.1:0.1);
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

subplot(1,4,2)
for i=1:2
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    chi2mat = (diff'*diff).*covi;
    semilogx(rsub_arr, diff./e,'-','DisplayName',...
        strcat('\chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on
end

h=legend('show','Location','southeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
ylim([-2,2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Delta/\sigma', 'fontsize',15);

subplot(1,4,3)
for i=3:4
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    semilogx(rsub_arr, diff./e,'-','DisplayName',...
        strcat('\chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on
end
h=legend('show','Location','southeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
ylim([-2,2])
ylabel('\Delta/\sigma', 'fontsize',15);

subplot(1,4,4)
for i=5:6
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    semilogx(rsub_arr, diff./e,'-','DisplayName',...
        strcat('\chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on
end
h=legend('show','Location','southeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
ylim([-2,2])
ylabel('\Delta/\sigma', 'fontsize',15);

savename = strcat(pltsavedir,'toy_chi2'); 
print(savename,'-dpng');close
%% toy model chi2 matrix
flight=40030;
inst=1;
mypaths=get_paths(flight);
pltsavedir = strcat(mypaths.alldat,'plots/cov_test/');
ifield = 4;
masklim = false;
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst));
load(sprintf('%s/stackdat_%s',loaddir,dt.name),'stackdatall');
load(sprintf('%s/excessdat_%s',loaddir,dt.name),'excessdatall');

figure
setwinsize(gcf,1500,600)
im = 2;

stackdat = stackdatall(im).stackdat;
excessdat = excessdatall(im).excessdat;
m_min = stackdat.m_min;
m_max = stackdat.m_max;
rsub_arr = stackdat.rsub_arr;

cov_mat = excessdat.excov.covcbsub;
covi = inv(cov_mat);
e = sqrt(diag(cov_mat))';

diff_list = [];
diff = e;
diff_list = [diff_list;diff];
diff = e;
diff(2:2:end) = -diff(2:2:end);
diff_list = [diff_list;diff];

for i=1:2
    diff = diff_list(i,:);
    chi2indept = sum(diff.^2./e.^2);
    chi2full = (diff*covi*diff');
    chi2mat = (diff'*diff).*covi;
    subplot(2,3,1)
    semilogx(rsub_arr, diff./e,'-');hold on
    
    subplot(2,3,2)
    semilogx(rsub_arr, diff,'-','DisplayName',strcat('case',num2str(i),...
        ': \chi^2_{indept} =',num2str(chi2indept,'%.2f'),...
        ', \chi^2_{full} =',num2str(chi2full,'%.2f')));hold on 
    
    subplot(2,3,3+i)
    imageclip(chi2mat);
    caxis([-7,7]);
    title(strcat('(Case ',num2str(i),') \Delta_i C^{-1}_{ij} \Delta_j'));
    xticks([3:4:15])
    xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    xtickangle(45)
    yticks([3:4:15])
    yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
        num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
    ytickangle(45)

    subplot(2,3,6)
    cum_chi2 = zeros(size(rsub_arr));
    for irsub=1:numel(rsub_arr)
        cum_chi2(irsub) = sum(sum(chi2mat(1:irsub,1:irsub)))/irsub;
    end
    loglog(rsub_arr, cum_chi2,'.-','DisplayName',strcat('case',num2str(i)));hold on
    
end

subplot(2,3,3)
imageclip(covi);
caxis([-1.5,1.5]);
h = colorbar;
ylabel(h, '[nW/m^2/sr]^{-2}');
title('C^{-1}');
xticks([3:4:15])
xticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
xtickangle(45)
yticks([3:4:15])
yticklabels({num2str(rsub_arr(3),'%.1e'),num2str(rsub_arr(7),'%.1e'),...
    num2str(rsub_arr(11),'%.1e'),num2str(rsub_arr(15),'%.1e')});
ytickangle(45)
if im==1
    ylabel('Excess Cov Inv.', 'fontsize',15);
end

subplot(2,3,1)
xlim([4e-1,1.1e3])
ylim([-1.2,1.2])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Delta/\sigma', 'fontsize',15);

subplot(2,3,2)
h=legend('show','Location','northeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Delta [nW/m^2/sr]', 'fontsize',12);

subplot(2,3,6)
h=legend('show','Location','northeast');
set(h,'fontsize',10)
xlim([4e-1,1.1e3])
xlabel('r [arcsec]', 'fontsize',15);
ylabel('\Sigma \Delta C^{-1} \Delta (<r) / DoF (<r)', 'fontsize',12);
savename = sprintf('%s%s_bgsubprofs',pltsavedir,dt.name);

savename = strcat(pltsavedir,'toy_cov'); 
print(savename,'-dpng');close
