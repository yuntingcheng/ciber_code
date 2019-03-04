%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%36277, 40030 read noise consistency check
% Use 40030 elat10 & 36277 BootesB, both has nfr=25
% compare their 1D and 2D PS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixsize=7;
iter_clip=3;
iter_mask=5;
%%
darkdir='/Users/ytcheng/ciber/doc/20150810_DarkCurrent/';
%savedir=strcat('/Users/ytcheng/ciber/doc/20160906_NoiseRealization/',...
%            'darkstat/',num2str(flight),'/');
        
savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_noise_flights/';
ell = get_l(1024,1024,pixsize,1);
%%
%%%%%%%%%%%%%%%%%%%%%%% halves 

dir1=sprintf('%s/40030/TM1/elat10/first/data/',darkdir);
dir2=sprintf('%s/40030/TM1/elat10/second/data/',darkdir);


%dir1=sprintf('%s/36277/TM1/BootesB/first/data/',darkdir);
%dir2=sprintf('%s/36277/TM1/BootesB/second/data/',darkdir);


scan1=dir(strcat(dir1,'dark*'));
scan2=dir(strcat(dir2,'dark*'));

Cld_arr=zeros(numel(scan1),29);
Cl2d_arr=zeros(numel(scan1),1024,1024);

for i=1:numel(scan1)
load(strcat(dir1,scan1(i).name));dark1=darkmap;
load(strcat(dir2,scan2(i).name));dark2=darkmap;

dark=(dark1-dark2)./sqrt(2);
[~,mdark]=get_skymap(dark,ones(1024),iter_clip);
dark=dc_offset_remove(dark,mdark);

[Cl,l,~,~,~,~,Cl2d]=get_angular_spec(dark.*mdark,dark.*mdark,pixsize);
Cld_arr(i,:)=Cl;
Cl2d_arr(i,:,:)=Cl2d;

%disp(sprintf('%d,TM%d,%s,%d/%d diff done',...
%                flight,inst,fieldname,i,numel(scanfull)));



figure
setwinsize(gcf,1500,300)

m1=subplot(1,3,1);
imageclip(dark.*mdark);

m2=subplot(1,3,2);
imagesc(log10(abs(Cl2d)));
colorbar
caxis([-12,-9])
%title(sprintf('SWIREdiff'));
colormap(m2,'default')

m3=subplot(1,3,3);
imagesc(log10(abs(Cl2d)));
xlim([513-50,513+50]); ylim([513-50,513+50]);
colorbar
caxis([-12,-9])
hold on
[c,h]=contour(ell,[500,1000,3000,5000,8000]);
clabel(c,h)
colormap(m3,'default')

end
%%
Cl2d4_arr=Cl2d_arr;
Cld4_arr=Cld_arr;

%Cl2d3_arr=Cl2d_arr;
%Cld3_arr=Cld_arr;
%%

figure
setwinsize(gcf,1500,300)

subplot(1,3,1);
imageclip(squeeze(log10(mean(Cl2d4_arr))));
caxis([-11,-9])
title('40030 mean 2DPS')

subplot(1,3,2);
imageclip(squeeze(log10(mean(Cl2d3_arr))));
caxis([-11,-9])
title('36277 mean 2DPS')

subplot(1,3,3);
imageclip(squeeze(mean(Cl2d4_arr))./squeeze(mean(Cl2d3_arr)));
title('40030/36277 mean 2DPS ratio')

print(strcat(savedir,'ave2DPS'),'-dpng');
%%
cp=get_cal_params('flight',36277);
cal=cp(1).apf2eps.*cp(1).eps2nWpm2ps;

plt3=loglog(l,l.*(l+1).*prctile(Cld3_arr.*cal.^2,50)./2./pi,'bo');hold on
errorbar(l,l.*(l+1).*prctile(Cld3_arr.*cal.^2,50)./2./pi,...
    l.*(l+1).*(prctile(Cld3_arr.*cal.^2,50)-prctile(Cld3_arr.*cal.^2,16))./2./pi,...
    l.*(l+1).*(prctile(Cld3_arr.*cal.^2,84)-prctile(Cld3_arr.*cal.^2,50))./2./pi,'bo');

plt4=loglog(l,l.*(l+1).*prctile(Cld4_arr.*cal.^2,50)./2./pi,'ro');hold on
errorbar(l,l.*(l+1).*prctile(Cld4_arr.*cal.^2,50)./2./pi,...
    l.*(l+1).*(prctile(Cld4_arr.*cal.^2,50)-prctile(Cld4_arr.*cal.^2,16))./2./pi,...
    l.*(l+1).*(prctile(Cld4_arr.*cal.^2,84)-prctile(Cld4_arr.*cal.^2,50))./2./pi,'ro');

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW m^{-2} sr^{-1})$',...
        'interpreter','latex','fontsize',18)
legend([plt3,plt4],{'26277','40030'},...
        'Location','southeast','FontSize',15);
legend boxoff
xlim([2e2,2e5]);
ylim([3e-2,1e4]);

print(strcat(savedir,'PS1D'),'-dpng');