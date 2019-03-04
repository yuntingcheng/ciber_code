%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get average signal Cl to put in simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
load(strcat(savedirfull,'flightCldat'),'flightCldat');
%% flight PS of all fields
numsig=zeros(1,21);
densig=zeros(1,21);
fig=figure;
for ifield=[8,7,6,5,4,2,1]
dt=get_dark_times(flight,inst,ifield);
l=flightCldat(ifield).l;
Cl_blcor=flightCldat(ifield).Cl_blcor;
Cl_blcor_err=flightCldat(ifield).Cl_blcor_err;
FFsigbias=flightCldat(ifield).FFsigbias;

%%%%% Raw flight + FF error correction factor + debias +bl correction%%%%%

yvalue=(l.*(l+1).*Cl_blcor./FFsigbias./2./pi);
stdtot=Cl_blcor_err./FFsigbias;
ebar=(l.*(l+1).*stdtot./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
errorbar(l,yvalue,ebarlow,ebar,'-o','color',get_color(ifield),...
    'markersize',3,'DisplayName',dt.name);hold on

numsig=numsig+Cl_blcor./FFsigbias./stdtot.^2;
densig=densig+1./stdtot.^2;
end
avgClsig=numsig./densig;

loglog(l,l.*(l+1).*avgClsig./2./pi,...
    'k','DisplayName','avg','linewidth',3);
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

save(strcat(savedirfull,'avgClsig'),'avgClsig');
