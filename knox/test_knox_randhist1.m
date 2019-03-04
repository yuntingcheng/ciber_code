loaddir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/rand_level/';
pixscale=7;
[Cl,~,~,l,~,dClknox,binl,dl] = Cl_from_Cl2d_clip(ones(1024),pixscale);


load(strcat(loaddir,'hist/data1'),'data');
load(strcat(loaddir,'hist/vpixuse_arr'),'vpixuse_arr');
%%
ell_arr=[20,26,29];
Cl2din=data(1).inavg_arr;

for isig=1:5
sig=data(isig).sig;
inavg_arr=data(isig).inavg_arr;
inavg2_arr=data(isig).inavg2_arr;
outavg_arr=data(isig).outavg_arr;
outavg2_arr=data(isig).outavg2_arr;
xpix_arr=data(isig).xpix_arr;
ypix_arr=data(isig).ypix_arr;
in_arr=data(isig).in_arr;
out_arr=data(isig).out_arr;
outstd=sqrt(outavg2_arr-outavg_arr.^2); 
for iell=1:3
    ell=ell_arr(iell);
    subxpix_arr=xpix_arr(iell*6-5:iell*6);
    subypix_arr=ypix_arr(iell*6-5:iell*6);
    subin_arr=in_arr(:,iell*6-5:iell*6);
    subout_arr=out_arr(:,iell*6-5:iell*6);
    subvin_arr=vpixuse_arr(iell*6-5:iell*6);
    figure
    setwinsize(gcf,1200,600)
    for ipix=1:6
    xpix=subxpix_arr(ipix);ypix=subypix_arr(ipix);
    
    subplot(2,3,ipix)
    h1=histogram(subout_arr(:,ipix));hold on
    h2=histogram(subin_arr(:,ipix),'BinEdges',h1.BinEdges);hold off
    legend([h2,h1],...
    {sprintf('input,avg=%.2e,std=%.2e',...
    mean(subin_arr(:,ipix))/Cl2din(xpix,ypix),std(subin_arr(:,ipix))/Cl2din(xpix,ypix)),...
    sprintf('output,avg=%.2e,std=%.2e',...
    mean(subout_arr(:,ipix))/Cl2din(xpix,ypix),std(subout_arr(:,ipix))/Cl2din(xpix,ypix))},...
       'location','northeast');
    if ipix==1
        title(sprintf('[%d,%d] max(Cl2Din)',xpix,ypix));
    elseif ipix==2
        title(sprintf('[%d,%d] med(Cl2Din)',xpix,ypix));
    elseif ipix==3
        title(sprintf('[%d,%d] min(Cl2Din)',xpix,ypix));   
    elseif ipix==4
        vout=outavg_arr(xpix,ypix)./Cl2din(xpix,ypix);
        vin=subvin_arr(ipix);
        title(sprintf('[%d,%d] max(<Cl2Dmap>/Cl2din)=%.2f',xpix,ypix,vin));
    elseif ipix==5
        vout=outavg_arr(xpix,ypix)./Cl2din(xpix,ypix);
        vin=subvin_arr(ipix);
        title(sprintf('[%d,%d] med(<Cl2Dmap>/Cl2din)=%.2f',xpix,ypix,vin));
    elseif ipix==6
        vout=outstd(xpix,ypix)./Cl2din(xpix,ypix);
        vin=subvin_arr(ipix);
        title(sprintf('[%d,%d] max(std(Cl2Dmap)/Cl2din)=%.2f',xpix,ypix,vin));
    end

    end
    savename=sprintf('%shist/ell%dsig%d',savedir,ell,isig);
    print(savename,'-dpng');close
end

end
