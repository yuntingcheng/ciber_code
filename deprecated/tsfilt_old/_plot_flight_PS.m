flight=40030;
inst=1;
savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
load(strcat(savedirfull,'flightCldat'),'flightCldat');
%% flight PS of individual field
for ifield=[8,7,6,5,4,2,1]
dt=get_dark_times(flight,inst,ifield);
l=flightCldat(ifield).l;
bl=flightCldat(ifield).bl;
Clflightf=flightCldat(ifield).Clflgihtf;
nClf_arr=flightCldat(ifield).nClf_arr;
FFcorfac=flightCldat(ifield).FFcorfac;
FFcorfacerr=flightCldat(ifield).FFcorfacerr;
Cl_blcor=flightCldat(ifield).Cl_blcor;
Cl_blcor_err=flightCldat(ifield).Cl_blcor_err;
fig=figure;
%%%%%% 'Raw flight' (no FF err correction, no de-bais) %%%%%%
yvalue=(l.*(l+1).*Clflightf./2./pi);
pltraw=plot(l,yvalue,'or','markersize',3);hold on

%%%%%% Raw flight debias %%%%%%
yvalue=(l.*(l+1).*(Clflightf-mean(nClf_arr))./2./pi);
ebar=(l.*(l+1).*std(nClf_arr)/2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltrawdb=errorbar(l,yvalue,ebarlow,ebar,'.r','markersize',10);hold on

%%%%%% Raw flight + FF error correction factor %%%%%%%%%%

yvalue=(l.*(l+1).*(Clflightf-FFcorfac)./2./pi);
ebar=(l.*(l+1).*FFcorfacerr./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltscale=errorbar(l,yvalue,ebarlow,ebar,'ob','markersize',3);hold on

%%%%%% Raw flight + FF error correction factor + debias%%%%%%%%%%

yvalue=(l.*(l+1).*(Clflightf-FFcorfac-mean(nClf_arr))./2./pi);
stdtot=sqrt(FFcorfacerr.^2+var(nClf_arr));
ebar=(l.*(l+1).*stdtot./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltscaledb=errorbar(l,yvalue,ebarlow,ebar,'.b','markersize',10);hold on

%%%%% Raw flight + FF error correction factor + debias +bl correction%%%%%

yvalue=(l.*(l+1).*(Cl_blcor)./2./pi);
stdtot=Cl_blcor_err;
ebar=(l.*(l+1).*stdtot./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltscaledb_bl=errorbar(l,yvalue,ebarlow,ebar,'om','markersize',3);hold on

%%%%%%%%%% RN +ph %%%%%%%%%%%%%%
y1=(l.*(l+1).*(prctile(nClf_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(nClf_arr,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on

xlim([1e2,2e5]);ylim([1e-2,1e4]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title(strcat(dt.name,'--full'));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
legend([pltraw,pltrawdb,pltscale,pltscaledb,pltscaledb_bl,pltn],...
 {'flight','flight-noise','flight-correction',...
 'flight-correction-noise','(flight-correction-noise)/bl','noise'},...
       'Location','southeast','FontSize',10);
legend boxoff
imname=strcat(savedirfull,'flightfull_i',num2str(ifield));
%print(imname,'-dpng');close
end

%% flight PS of all fields
fig=figure;
for ifield=[8,7,6,5,4,2,1]
dt=get_dark_times(flight,inst,ifield);
l=flightCldat(ifield).l;
bl=flightCldat(ifield).bl;
Clflightf=flightCldat(ifield).Clflgihtf;
nClf_arr=flightCldat(ifield).nClf_arr;
FFcorfac=flightCldat(ifield).FFcorfac;
FFcorfacerr=flightCldat(ifield).FFcorfacerr;
Cl_blcor=flightCldat(ifield).Cl_blcor;
Cl_blcor_err=flightCldat(ifield).Cl_blcor_err;

%%%%% Raw flight + FF error correction factor + debias +bl correction%%%%%

yvalue=(l.*(l+1).*Cl_blcor./2./pi);
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
imname=strcat(savedirfull,'allfields');
%print(imname,'-dpng');close
