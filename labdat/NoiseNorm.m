flight=40030;
inst=1;
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
g1=cp(inst).apf2eps;g2=cal./g1;

savedir=strcat('/Users/ytcheng/ciber/doc/20161205_NoiseNorm/TM',...
        num2str(inst),'/');
halfdatdir='/Users/ytcheng/ciber/doc/20160920_FlightDiff/';
load(strcat(halfdatdir,'band',num2str(inst),'_halfdat'),'halfdat');
load(strcat(halfdatdir,'band',num2str(inst),'_diffpsdat'),'psdat');
%% get the field and sim PS
%{
for ifield=4:8
    mask=halfdat(ifield).bigmask;
    mkk=halfdat(ifield).wMkk;
    fielddiff=(halfdat(ifield).raw1-halfdat(ifield).raw2)./sqrt(2);
    [Cld,~,Cl2dd,l]=get_Cl(fielddiff,mask,mkk,pixscale,ones(1024));
    nsimps=noise_sim_ps(flight,inst,ifield,mask,mkk,ones(1024),1,g1);
    data(ifield).nsimps=nsimps;
    data(ifield).field.Cl2dd=Cl2dd;
    data(ifield).field.Cld=Cld;
end
save(strcat(savedir,'TM',num2str(inst),'_data'),'data');
%}
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate chi2
load(strcat(savedir,'TM',num2str(inst),'_data'),'data');
[~,~,~,l]=get_Cl(randn(1024),ones(1024),eye(29),pixscale,ones(1024));
Arn_arr=0.2:0.01:1;
chi2_arr=zeros(1,numel(Arn_arr));
for iamp=1:numel(Arn_arr)
Arn=Arn_arr(iamp);
nave_arr=zeros(16,5);
nstd_arr=zeros(16,5);
field_arr=zeros(16,5);
ell_arr=zeros(16,5);
for i=4:8
rnCl_arr=data(i).nsimps.rnCl_arr;
if inst==1
phCl_arr=data(i).nsimps.phCl_arr.*(g1/-2.68);
else
phCl_arr=data(i).nsimps.phCl_arr.*(g1/-2.88);
end
nCl_arr=(phCl_arr+rnCl_arr.*Arn).*cal^2;
nave=mean(nCl_arr);nave([1,2,12,end-1,end])=[];
nstd=std(nCl_arr);nstd([1,2,12,end-1,end])=[];
nave_arr(:,i-3)=nave;
nstd_arr(:,i-3)=nstd;
field=data(i).field.Cld;field([1,2,12,end-1,end])=[];
field_arr(:,i-3)=field;
ell=l;ell([1,2,12,end-1,end])=[];
ell_arr(:,i-3)=ell;
end
nave_arr=nave_arr(:)';
nstd_arr=nstd_arr(:)';
field_arr=field_arr(:)';
ell_arr=ell_arr(:)';
chi2=sum((field_arr-nave_arr).^2./(nstd_arr).^2);

chi2_arr(iamp)=chi2;


if iamp==46
figure
loglog(ell_arr,((field_arr-nave_arr)./nstd_arr).^2,'ro');
hold on
end

if Arn==1
loglog(ell_arr,((field_arr-nave_arr)./nstd_arr).^2,'bo');
xlim([1e2,2e5]);ylim([1e-5,1e3]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\chi^2$','interpreter','latex','fontsize',18)
legend({'Arn=0.65','Arn=1'},'location','southeast')
print(strcat(savedir,'chi2_ell'),'-dpng');%close
end

end
%% plot Arn - chi2/dof relation
chi2min=chi2_arr(find(chi2_arr==min(chi2_arr)));
Arnmin=Arn_arr(find(chi2_arr==min(chi2_arr)));
ndof=numel(field_arr)-1;
semilogy(Arn_arr,chi2_arr./ndof,'.');
xlabel('A_{RN}')
ylabel('$\chi^2/dof$','interpreter','latex','fontsize',18)
title(sprintf('min(Chi2)=%.3e,Chi2/dof=%.2f,Arn=%.2f',...
                            chi2min,chi2min/ndof,Arnmin))
print(strcat(savedir,'chi2_Arn'),'-dpng');%close
%%
for i=4:8
rnCl_arr=data(i).nsimps.rnCl_arr;
phCl_arr=data(i).nsimps.phCl_arr;

nClorig_arr=(phCl_arr+rnCl_arr).*cal^2;
if inst==1
nCl_arr=(phCl_arr.*(g1/-2.68)+rnCl_arr.*Arnmin).*cal^2;
else
nCl_arr=(phCl_arr.*(g1/-2.88)+rnCl_arr.*Arnmin).*cal^2;
end

Cld=data(i).field.Cld;

fig=figure;
errorbar(l,l.*(l+1).*prctile(nCl_arr,50)./2./pi,...
     l.*(l+1).*(prctile(nCl_arr,50)-prctile(nCl_arr,16))./2./pi,...
     l.*(l+1).*(prctile(nCl_arr,84)-prctile(nCl_arr,50))./2./pi,...
                    '.','color','b','markersize',10);hold on
                
errorbar(l,l.*(l+1).*prctile(nClorig_arr,50)./2./pi,...
     l.*(l+1).*(prctile(nClorig_arr,50)-prctile(nClorig_arr,16))./2./pi,...
     l.*(l+1).*(prctile(nClorig_arr,84)-prctile(nClorig_arr,50))./2./pi,...
                    '.','color',[0.6,0.6,0.6],'markersize',10);hold on

plot(l,l.*(l+1).*Cld./2./pi,'r+','markersize',10);

xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
title(halfdat(i).name)
xlim([1e2,2e5]);

if inst==1
ylim([1e-2,1e4]);
else
ylim([5e-3,2e3]);
end

ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log')
legend({sprintf('Arn=%.2f',Arnmin),'Arn=1','flight diff'},...
       'Location','southeast','FontSize',15);
legend boxoff
print(strcat(savedir,'PSi',num2str(i)),'-dpng');%close
end
