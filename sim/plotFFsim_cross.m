%%
savedir='/Users/ytcheng/ciber/doc/20170320_FFsim/xCl_sim/';
loaddir='/Users/ytcheng/ciber/doc/20170320_FFsim/';
load(strcat(loaddir,'simCldat'),'simCldat');
load(strcat(loaddir,'simsigCldat1'),'simsigCldat');

l=simsigCldat.l;
for i=1:21
    fig=figure;
    Clx=simsigCldat.cross(i).meanClx;
    Clxerr=simsigCldat.cross(i).stdClx;
    
    Clx1=Clx-simCldat.cross(i).meanClx;
    Clx1err=sqrt(simCldat.cross(i).stdClx.^2+Clxerr.^2);

    Clxpos=Clx(find(Clx>0));lpos=l(find(Clx>0));
    Clxneg=Clx(find(Clx<0));lneg=l(find(Clx<0));
    Clxerrpos=Clxerr(find(Clx>0));
    Clxerrneg=Clxerr(find(Clx<0));
    
    Clx1pos=Clx1(find(Clx1>0));l1pos=l(find(Clx1>0));
    Clx1neg=Clx1(find(Clx1<0));l1neg=l(find(Clx1<0));
    Clx1errpos=Clx1err(find(Clx1>0));
    Clx1errneg=Clx1err(find(Clx1<0));
    
    % slightly offset for better plot presentation
    l1pos=l1pos.*1.05;l1neg=l1neg.*1.05;

    errorbar(lneg,lneg.*(lneg+1).*-Clxneg./2./pi,...
        lneg.*(lneg+1).*Clxerrneg./2./pi,'b.','markersize',10);hold on    
    errorbar(lpos,lpos.*(lpos+1).*Clxpos./2./pi,...
        lpos.*(lpos+1).*Clxerrpos./2./pi,'r.','markersize',10);hold on

    errorbar(l1neg,l1neg.*(l1neg+1).*-Clx1neg./2./pi,...
        l1neg.*(l1neg+1).*Clx1errneg./2./pi,'bo');    
    errorbar(l1pos,l1pos.*(l1pos+1).*Clx1pos./2./pi,...
        l1pos.*(l1pos+1).*Clx1errpos./2./pi,'ro');
    
    ax = get(fig,'CurrentAxes');
    set(ax,'XScale','log','YScale','log');
   
    name1=simsigCldat.siminfo(simsigCldat.cross(i).field(1)).name;
    name2=simsigCldat.siminfo(simsigCldat.cross(i).field(2)).name;
    title(sprintf('%s x %s',name1,name2));
    xlim([[1e2,2e5]]);
    ylim([1e-2,1e3]);
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\ell(\ell+1)|C_\ell|/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
    legend({'x-full sim(-)','x-full sim(+)',...
'x-full sim - x-noise(-)','x-full sim - x-noise(+)'},...
        'location','southeast','fontsize',10)
    legend boxoff
    drawnow
    imname=strcat(savedir,'xCl',num2str(i));
    print(imname,'-dpng');close

end

