%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adhoc read noise model validation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=36277;
inst=1;
pixscale=7;
iter_mask=10;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
%%
DCdir=strcat('/Users/ytcheng/ciber/doc/20160808_DarkProcess/',...
            num2str(flight),'/');
loaddir=strcat('/Users/ytcheng/ciber/doc/20160906_NoiseRealization/',...
            'darkstat/',num2str(flight),'/');
load(strcat(loaddir,'TM',num2str(inst),'_darkps'),'darkps');
savedir=strcat('/Users/ytcheng/ciber/doc/20160906_NoiseRealization/',...
    'validation/',num2str(flight),'/');
%% 
type='diff';
for i=1:numel(darkps)

Clsim_arr=zeros(100,29);
for j=1:100   
if type=='full'    
    rnmap=readnoise_realization(darkps(i).Clf2d_ave,pixscale);
else
    rnmap=readnoise_realization(darkps(i).Cld2d_ave,pixscale);
end

Clsim=get_angular_spec(rnmap,rnmap,pixscale);
Clsim_arr(j,:)=Clsim;
end

if type=='full'
    Cl_arr=darkps(i).Clf_arr;
else
    Cl_arr=darkps(i).Cld_arr;
end

l=darkps(i).l;
Clsim_arr=Clsim_arr.*cal.^2;Cl_arr=Cl_arr.*cal.^2;

%%% plot PS 
figure
pltsim=loglog(l,l.*(l+1).*prctile(Clsim_arr,50)./2./pi,'ro');hold on
errorbar(l,l.*(l+1).*prctile(Clsim_arr,50)./2./pi,...
    l.*(l+1).*(prctile(Clsim_arr,50)-prctile(Clsim_arr,16))./2./pi,...
    l.*(l+1).*(prctile(Clsim_arr,84)-prctile(Clsim_arr,50))./2./pi,'ro');

pltdat=loglog(l,l.*(l+1).*prctile(Cl_arr,50)./2./pi,'bo');hold on
errorbar(l,l.*(l+1).*prctile(Cl_arr,50)./2./pi,...
    l.*(l+1).*(prctile(Cl_arr,50)-prctile(Cl_arr,16))./2./pi,...
    l.*(l+1).*(prctile(Cl_arr,84)-prctile(Cl_arr,50))./2./pi,'bo');
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW m^{-2} sr^{-1})$',...
        'interpreter','latex','fontsize',18)
legend([pltdat,pltsim],{'real dark data','simulated dark data'},...
        'Location','southeast','FontSize',15);
legend boxoff
title(sprintf('Noise Test for TM%d, %d, %s',...
                 inst,flight,darkps(i).field));    
xlim([2e2,2e5]);

if type=='full'
%ylim([1e-2,5e3]);%TM1
ylim([1e-3,5e2]);%TM2
else 
%ylim([1e-1,1e4]);%TM1
ylim([1e-2,1e3]);%TM2
end

drawnow

savename=sprintf('%sTM%d_i%i_PS_%s',savedir,inst,i,type);
print(savename,'-dpng');%close

%%% plot PS ratio
norm=prctile(Cl_arr,50);
figure
pltsim=semilogx(l,prctile(Clsim_arr,50)./norm,'ro');hold on
errorbar(l,prctile(Clsim_arr,50)./norm,...
    (prctile(Clsim_arr,50)-prctile(Clsim_arr,16))./norm,...
    (prctile(Clsim_arr,84)-prctile(Clsim_arr,50))./norm,'ro');
pltdat=semilogx(l,prctile(Cl_arr,50)./norm,'bo');hold on
errorbar(l,prctile(Cl_arr,50)./norm,...
    (prctile(Cl_arr,50)-prctile(Cl_arr,16))./norm,...
    (prctile(Cl_arr,84)-prctile(Cl_arr,50))./norm,'bo');
xlim([2e2,2e5]);ylim([-3,7]);
plot(xlim,[1,1],':k')

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$C_\ell/C_\ell^{data}$',...
        'interpreter','latex','fontsize',18)
legend([pltdat,pltsim],{'real dark data','simulated dark data'},...
        'Location','southeast','FontSize',15);
legend boxoff
title(sprintf('Noise Model Ratio for TM%d, %d, %s',...
                    inst,flight,darkps(i).field));    
drawnow

savename=sprintf('%sTM%d_i%i_PS_%s_ratio',savedir,inst,i,type);
print(savename,'-dpng');%close
end
                 
