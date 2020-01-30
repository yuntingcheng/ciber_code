
figure
for ifield=4:8
m_arr=10:0.5:23;
rad_arr = get_mask_radius_th(1,ifield,m_arr,1);
plot(m_arr,rad_arr);hold on
end
legend({'4','5','6','7','8'})

%%
flight = 40030;
mypaths=get_paths(flight);
m_min_arr = [4,4,4,4,12,13,15,16:19];
m_max_arr = [9,10,11,12,13,14,16,17:20];
Njack = 50;

figure
setwinsize(gcf, 800, 600)

for inst=1:2
bandname = ['I','H'];
pltsavedir=(strcat(mypaths.alldat,'plots/TM',num2str(inst),'/'));
names={};
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

for ifield = 4:8

    dt=get_dark_times(flight,inst,ifield);
    names{end+1}=dt.name;
    load(sprintf('%s/psfdat_%s',savedir,dt.name),'psfdatall');

    for im=1
        psfdat = psfdatall.comb(im);
        r_arr = psfdat.r_arr;
        p = psfdat.all.profcb;
        im_in = psfdat.im_in;
        im_mid = psfdat.im_mid;
        im_out = psfdat.im_out;
        dat_profcb = zeros([Njack,numel(r_arr)]);
        for isub=1:Njack
            dat_profcb(isub,:) = psfdat.jack(isub).profcb;
        end
        covcb = get_cov_matrix(dat_profcb).*(Njack-1);
        e = sqrt(diag(covcb)');
        
        subplot(2,2,2*inst-1)
        p1=semilogx(r_arr, p,'o-','linewidth',2,'markersize',8,...
        'color',get_color(ifield-3));hold on
        subplot(2,2,2*inst)
        loglog(r_arr, p,'o-','linewidth',2,'markersize',8,...
        'color',get_color(ifield-3));hold on

    end
end

subplot(2,2,2*inst-1)
ylim([3e-7,1.02])
xlim([4e-1,1e2])
p=vline(7,'k:');
set(p,'linewidth',3);
grid on
xlabel('r [arcsec]', 'fontsize',20)
ylabel('PSF', 'fontsize',20)
legend(names,'fontsize',12);
subplot(2,2,2*inst)
ylim([3e-5,1.1])
xlim([4e-1,1e2])
p=vline(7,'k:');
set(p,'linewidth',3);
grid on
xlabel('r [arcsec]', 'fontsize',20)
ylabel('PSF', 'fontsize',20)
end
savename = sprintf('%s%s/psf',mypaths.alldat,'plots/paperfig/');
print(savename,'-dpng');close