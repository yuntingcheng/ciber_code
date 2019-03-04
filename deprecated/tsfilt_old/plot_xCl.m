flight=40030;
inst=1;
pixscale=7;
savedir='/Users/ytcheng/ciber/doc/20170325_alldat/TM1/xCl_flight/';
load('/Users/ytcheng/ciber/doc/20170325_alldat/TM1/xmkk','xmkk');
load('/Users/ytcheng/ciber/doc/20170320_FFsim/simCldat','simCldat');
load('/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/fwfulldat',...
                                                             'fwfulldat');
for ifield=[8,7,6,5,4,2,1]
dt=get_dark_times(flight,inst,ifield);
loaddir=strcat('/Volumes/HD1TB/CIBER/tsfilt/',...
       'sinfiltamp/FullMap/field',num2str(ifield),'/');
load(strcat(loaddir,'flightmap'),'flightmap');
load(strcat(loaddir,'maskin'),'maskin');
std_noise=fwfulldat(ifield).std_fCl2d;
weight=(fftshift(fftshift(1./std_noise)))';

map=flightmap.calmapf;
map=(map-mean(map(find(map)))).*maskin;
xinfo(ifield).name=dt.name;
xinfo(ifield).mask=maskin;
xinfo(ifield).weight=weight;
xinfo(ifield).map=map;
end
%%
count=0;
for ifield=[8,7,6,5,4,2,1]
    maski=xinfo(ifield).mask;
    weighti=xinfo(ifield).weight;
    mapi=xinfo(ifield).map;
    for jfield=[8,7,6,5,4,2,1]
        if ifield>jfield
            count=count+1;
            maskj=xinfo(jfield).mask;
            weightj=xinfo(jfield).weight;
            mapj=xinfo(jfield).map;

            mask=maski.*maskj;
            weight=(weighti+weightj)./2;
            Mkk=xmkk(count).Mkk;
            [Clx,~,~,l]=get_Clx(mapi,mapj,mask,Mkk,pixscale,weight);
            Clx1=Clx-simCldat.cross(count).meanClx;
            Clx1err=simCldat.cross(count).stdClx;
            
            Clxpos=Clx(find(Clx>0));lpos=l(find(Clx>0));
            Clxneg=Clx(find(Clx<0));lneg=l(find(Clx<0));
            Clx1pos=Clx1(find(Clx1>0));l1pos=l(find(Clx1>0));
            Clx1neg=Clx1(find(Clx1<0));l1neg=l(find(Clx1<0));
            Clx1errpos=Clx1err(find(Clx1>0));
            Clx1errneg=Clx1err(find(Clx1<0));
            
            % slightly offset for better plot presentation
            l1pos=l1pos.*1.05;l1neg=l1neg.*1.05;

            figure
            loglog(lneg,lneg.*(lneg+1).*-Clxneg./2./pi,'b.',...
                'markersize',10);hold on            
            loglog(lpos,lpos.*(lpos+1).*Clxpos./2./pi,'r.',...
                'markersize',10);hold on
            
            errorbar(l1neg,l1neg.*(l1neg+1).*-Clx1neg./2./pi,...
                l1neg.*(l1neg+1).*Clx1errneg./2./pi,'bo');
            errorbar(l1pos,l1pos.*(l1pos+1).*Clx1pos./2./pi,...
                l1pos.*(l1pos+1).*Clx1errpos./2./pi,'ro');
            
            title(strcat(xinfo(ifield).name,' x ',xinfo(jfield).name));
            xlim([[1e2,2e5]]);
            ylim([1e-2,1e3]);
            xlabel('$\ell$','interpreter','latex','fontsize',18)
            ylabel('$\ell(\ell+1)|C_\ell|/2\pi(nW^2 m^{-4} sr^{-2})$',...
                'interpreter','latex','fontsize',18)
            legend({'cross flight(-)','cross flight(+)',...
      'cross flight - cross noise(-)','cross flight - cross noise(+)'},...
                'location','southeast','fontsize',10)
            legend boxoff
            drawnow
            imname=strcat(savedir,'xCl',num2str(count));
            print(imname,'-dpng');close
        end
    end
end
