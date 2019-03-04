%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%In some DC map, the bottom-left quadrant has some DC offset.
%This test compares 2 dark map: one has obvious DC offset,
%while the other doesn't. 
%These maps are arbitray chosen by looking at the DC map.
%The DC correction is applied in both case, and compare the output PS 
%before and after correction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

savedir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_dc_offset/';
indir='/Users/ytcheng/ciber/doc/20150810_DarkCurrent/TM2/DGL/full/data/';
load(strcat(indir,'dark_14-45-38_68_106.mat'));dark1=darkmap;
load(strcat(indir,'dark_14-23-09_3_41.mat'));dark2=darkmap;
load(strcat(indir,'dark_15-10-18_13_51.mat'));dark3=darkmap;

diff1=(dark1-dark3)./sqrt(2);
diff2=(dark2-dark3)./sqrt(2);

[~,mask1]=get_skymap(diff1,ones(1024),3);
[~,mask2]=get_skymap(diff2,ones(1024),3);

[~,l,~,~,lbin]=get_angular_spec(dark1,dark1,7);
mkk1=get_mkk_sim(mask1,7,lbin,10,numel(lbin),1,ones(1024),0,NaN);
mkk2=get_mkk_sim(mask2,7,lbin,10,numel(lbin),1,ones(1024),0,NaN);

[cCl1,wCl1,Cl2d1,l,binl]=get_Cl(diff1,mask1,mkk1,7,ones(1024));
[cCl2,wCl2,Cl2d2,l,binl]=get_Cl(diff2,mask2,mkk2,7,ones(1024));
diff11=dc_offset_remove(diff1,mask1);
diff22=dc_offset_remove(diff2,mask2);
[cCl11,wCl11,Cl2d11,l,binl]=get_Cl(diff11,mask1,mkk1,7,ones(1024));
[cCl22,wCl22,Cl2d22,l,binl]=get_Cl(diff22,mask2,mkk1,7,ones(1024));
%% plot PS
cl1=loglog(l,l.*(l+1).*cCl1./2./pi,'b--');hold on
cl11=loglog(l,l.*(l+1).*cCl11./2./pi,'b');hold on
cl2=loglog(l,l.*(l+1).*cCl2./2./pi,'r--');hold on
cl22=loglog(l,l.*(l+1).*cCl22./2./pi,'r');hold on

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)

legend([cl1,cl11,cl2,cl22],...
{'dark1','dark1-DC shift','dark2','dark2-DC shift'},...
       'Location','southeast','FontSize',15);
legend boxoff
imname=strcat(savedir,'test_dc_offset_PS');
print(imname,'-dpng');%close
%% plot map
figure
subplot(2,2,1)
imageclip(diff1);
title('dark1');
subplot(2,2,2)
imageclip(diff11);
title('dark1-DC shift');
subplot(2,2,3)
imageclip(diff2);
title('dark2');
subplot(2,2,4)
imageclip(diff22);
title('dark2-DC shift');
imname=strcat(savedir,'test_dc_offset_map');
print(imname,'-dpng');%close
