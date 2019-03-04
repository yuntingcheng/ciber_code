function plot_flight_auto(flight,inst)
mypaths=get_paths(flight);
loaddir=strcat(mypaths.ciberdir,'doc/20170325_alldat/TM',num2str(inst),'/');
load(strcat(loaddir,'flightCldat'),'flightCldat');

savedir=strcat(mypaths.ciberdir,'doc/20170325_alldat/TM',...
    num2str(inst),'/plot_auto/');
%% flight PS of individual field
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    l=flightCldat(ifield).l;
    Clflight=flightCldat(ifield).Clflgiht;
    nCl_arr=flightCldat(ifield).nCl_arr;
    FFcorfac=flightCldat(ifield).FFcorfac;
    FFcorfacerr=flightCldat(ifield).FFcorfacerr;
    Cl_blcor=flightCldat(ifield).Cl_blcor;
    Cl_blcor_err=flightCldat(ifield).Cl_blcor_err;

    fig=figure;
    
    %%%%%% 'Raw flight' (no any debias) %%%%%%
    yvalue=(l.*(l+1).*Clflight./2./pi);
    plot(l,yvalue,'ok','markersize',5);hold on

    %%%%%%%%%% noise(RN +ph) %%%%%%%%%%%%%%
    y1=(l.*(l+1).*(prctile(nCl_arr,16))./2./pi);
    y2=(l.*(l+1).*(prctile(nCl_arr,84))./2./pi);
    pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
        [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on
        
%     %%%%%% Raw flight - FFbias %%%%%%%%%%
%     yvalue=(l.*(l+1).*(Clflight-FFcorfac)./2./pi);
%     ebar=(l.*(l+1).*FFcorfacerr./2./pi);
%     ebarlow=ebar;
%     for ii=1:numel(l)
%     if yvalue(ii)<=ebar(ii)
%         ebarlow(ii)=yvalue(ii)-1e-8;
%     end
%     end
%     errorbar(l,yvalue,ebarlow,ebar,'ob','markersize',5);hold on

    %%%%%% Raw flight -FFbias - noise %%%%%%%%%%
    yvalue=(l.*(l+1).*(Clflight-FFcorfac-mean(nCl_arr))./2./pi);
    stdtot=sqrt(FFcorfacerr.^2+var(nCl_arr));
    ebar=(l.*(l+1).*stdtot./2./pi);
    ebarlow=ebar;
    for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii)
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    end
    errorbar(l,yvalue,ebarlow,ebar,'ob','markersize',5);hold on

    %%%%% (Raw flight -FFbias - noise)/bl%%%%%
    yvalue=(l.*(l+1).*(Cl_blcor)./2./pi);
    stdtot=Cl_blcor_err;
    ebar=(l.*(l+1).*stdtot./2./pi);
    ebarlow=ebar;
    for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii)
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    end
    errorbar(l,yvalue,ebarlow,ebar,'.r','markersize',15);hold on
    
    %%%% shot noise ref line, fit the last 5 ell bins %%%%%
    const = mean(Cl_blcor(Cl_blcor>0 & l>3e4));
    plot(l,l.*(l+1).*const./2./pi,'r:','linewidth',1.5);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xlim([1e2,2e5]);ylim([1e-2,1e5]);
    ax = get(fig,'CurrentAxes');
    set(ax,'XScale','log','YScale','log');

    title(strcat(dt.name));
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
            'interpreter','latex','fontsize',18)
    h=legend({'flight','noise','flight-FFbias-noise',...
        '(flight-FFbias-noise)/bl','Cl=const fit'},'Location','southeast');
    set(h,'fontsize',10)
    legend boxoff

    imname=strcat(savedir,'flight_i',num2str(ifield));
    %print(imname,'-dpng');close
end
%% flight PS of all fields

fig=figure;
for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    l=flightCldat(ifield).l;
    bl=flightCldat(ifield).bl;
    Clflight=flightCldat(ifield).Clflgiht;
    nCl_arr=flightCldat(ifield).nCl_arr;
    FFcorfac=flightCldat(ifield).FFcorfac;
    FFcorfacerr=flightCldat(ifield).FFcorfacerr;
    Cl_blcor=flightCldat(ifield).Cl_blcor;
    Cl_blcor_err=flightCldat(ifield).Cl_blcor_err;


    %%%%% (Raw flight -FFbias - noise)/bl%%%%%
    yvalue=(l.*(l+1).*(Cl_blcor)./2./pi);
    stdtot=Cl_blcor_err;
    ebar=(l.*(l+1).*stdtot./2./pi);
    ebarlow=ebar;
    for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii)
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    end
    errorbar(l,yvalue,ebarlow,ebar,'-o','color',get_color(ifield),...
        'markersize',3,'DisplayName',dt.name);hold on
end

xlim([1e2,2e5]);ylim([1e-2,5e4]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title('all fields');
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','southeast');
set(h,'fontsize',10)
legend boxoff
imname=strcat(savedir,'flight_all');
print(imname,'-dpng');close

end