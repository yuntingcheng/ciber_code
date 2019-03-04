%%
savedir='/Users/ytcheng/ciber/doc/20170320_FFsim/';

load(strcat(savedir,'simCldat'),'simCldat');
load(strcat(savedir,'simsigCldat1'),'simsigCldat');
savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
load(strcat(savedirfull,'flightCldat'),'flightCldat');
%%
l=simsigCldat.l;
for ifield =[8,7,6,5,4,2,1]

fig=figure;

Clflightf=flightCldat(ifield).Clflgihtf;
yvalue=(l.*(l+1).*Clflightf./2./pi);
plot(l,yvalue,'.r','markersize',10,'DisplayName','flight');hold on

Cls=simsigCldat.auto(ifield).Clsc_arr;                
errorbar(l,l.*(l+1).*prctile(Cls,50)./2./pi,...
     l.*(l+1).*(prctile(Cls,50)-prctile(Cls,16))./2./pi,...
     l.*(l+1).*(prctile(Cls,84)-prctile(Cls,50))./2./pi,...
     '.k','markersize',10,'DisplayName','sim input sig * PSF');hold on

Cl1=simsigCldat.auto(ifield).Cl1_arr;                
errorbar(l,l.*(l+1).*prctile(Cl1,50)./2./pi,...
     l.*(l+1).*(prctile(Cl1,50)-prctile(Cl1,16))./2./pi,...
     l.*(l+1).*(prctile(Cl1,84)-prctile(Cl1,50))./2./pi,...
     'or','markersize',3,'DisplayName','sim output');hold on

Cln=flightCldat(ifield).nClf_arr;                
y1=(l.*(l+1).*(prctile(Cln,16))./2./pi);
y2=(l.*(l+1).*(prctile(Cln,84))./2./pi);
pltCln=fill([l,flip(l)],[abs(y1),abs(flip(y2))],[0.2,0.6,0.2],...
    'facealpha',0.5,'EdgeColor','none','DisplayName','RN+ph');hold on


Cl1dFF=mean(Cl1)-simCldat.auto(ifield).meanCld;
Cl1dFFerr=sqrt(simCldat.auto(ifield).stdCld.^2+var(Cl1));
yvalue=(l.*(l+1).*Cl1dFF./2./pi);
ebar=(l.*(l+1).*Cl1dFFerr./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
errorbar(l,yvalue,ebarlow,ebar,'ob','markersize',3,...
    'DisplayName','sig output - noiseFFbias');hold on

Cl1dFFn=Cl1dFF-mean(Cln);
Cl1dFFnerr=sqrt(Cl1dFFerr.^2+var(Cln));
yvalue=(l.*(l+1).*Cl1dFFn./2./pi);
ebar=(l.*(l+1).*Cl1dFFnerr./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
if yvalue(ii)<=ebar(ii)
    ebarlow(ii)=yvalue(ii)-1e-8;
end
end
pltdFF=errorbar(l,yvalue,ebarlow,ebar,'.b','markersize',10,...
    'DisplayName','sig output - noiseFFbias - noise');hold on

xlim([2e2,2e5]);ylim([1e-2,1e4]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
title(simsigCldat.siminfo(ifield).name); 
h=legend('show','Location','southeast');
set(h,'fontsize',10)
legend boxoff
savename=strcat(savedir,'FFsim_i',num2str(ifield));
print(savename,'-dpng');close    
end

