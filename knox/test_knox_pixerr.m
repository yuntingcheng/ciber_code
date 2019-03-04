loaddir='/Users/ytcheng/ciber/doc/20160808_DarkProcess/test_knox/rand_level/';
pixscale=7;
[~,~,~,l,~,dClknox,binl,dl,~,~,nell] = Cl_from_Cl2d_clip(ones(1024),pixscale);


load(strcat(loaddir,'hist/data1'),'data');data1=data;
load(strcat(loaddir,'hist/dataw'),'data');dataw=data;
load(strcat(loaddir,'hist/data'),'data');% for weight only
%%
for isig=1%1:5
sig=data(isig).sig;

weight=1./sqrt(data(isig).outavg2_arr-data(isig).outavg_arr.^2);
%weight=ones(1024);
dClrn=zeros(1,29);dClwn=zeros(1,29);

ell=get_l(1024,1024,pixscale,1);

for iell=9:29
mask=zeros(1024);mask((ell >= binl(iell)) & (ell <= binl(iell+1)))=1;
w_arr=weight(find(mask));

%%%% RN %%%%%%%%%
avg2_arr=data1(isig).outavg2_arr(find(mask));
avg_arr=data1(isig).outavg_arr(find(mask));
var_arr=avg2_arr-avg_arr.^2;
std_arr=sqrt(avg2_arr-avg_arr.^2);

meanCl=sum(w_arr.*avg_arr)./sum(w_arr);
varCl=sum(w_arr.^2.*var_arr)./sum(w_arr)^2;

dClrn(iell)=sqrt(varCl)/meanCl;
%%%% white %%%%%%%
avg2_arr=dataw(isig).outavg2_arr(find(mask));
avg_arr=dataw(isig).outavg_arr(find(mask));
var_arr=avg2_arr-avg_arr.^2;
std_arr=sqrt(avg2_arr-avg_arr.^2);

meanCl=sum(w_arr.*avg_arr)./sum(w_arr);
varCl=sum(w_arr.^2.*var_arr)./sum(w_arr)^2;

dClwn(iell)=sqrt(varCl)/meanCl;

end
figure
loglog(l,dClknox,'k');hold on
loglog(l,1./sqrt(nell),'k--');hold on
loglog(l,dClrn,'r.','markersize',10);
loglog(l,dClwn,'b.','markersize',10);
title(sprintf('sigma=%.2f',sig))
legend({'Knox','1/sqrt(N)','RN','white'})
xlabel('$\ell$','interpreter','latex','fontsize',18);
ylabel('$\delta C_\ell/\bar{C_\ell}$',...
        'interpreter','latex','fontsize',18);
end