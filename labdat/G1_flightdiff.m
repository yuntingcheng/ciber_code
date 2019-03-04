%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Cal G1 with flight diff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=2;
pixscale=7;
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

savedir='/Users/ytcheng/ciber/doc/20160912_CalFac/flightdiff/';
halfdatdir='/Users/ytcheng/ciber/doc/20160920_FlightDiff/';
load(strcat(halfdatdir,'band',num2str(inst),'_halfdat'),'halfdat');
alldatdir='/Users/ytcheng/ciber/doc/20160919_FlightDat/';
load(strcat(alldatdir,'band',num2str(inst),'_alldat'),'alldat');
%%
g1_arr=3:0.5:7;

for ig1=1:numel(g1_arr)
g1=-g1_arr(ig1);
chi2_arr=zeros(5,21);    
for ifield=4:8
mask=halfdat(ifield).bigmask;mkk=halfdat(ifield).wMkk;
rawdmap=halfdat(ifield).rawd;
nsimps=noise_sim_ps(flight,inst,ifield,mask,mkk,1,g1);
[psflight,~,cl2dflight]=get_Cl(rawdmap,mask,mkk,pixscale,nsimps.weight);

nmean=mean(nsimps.nCl_arr.*cal.*cal);nstd=std(nsimps.nCl_arr.*cal.*cal);
chi2=((psflight-nmean).^2)./(nstd.^2);
chi2_arr(ifield-3,:)=chi2;
chi2tot=sum(chi2(end-4:end));
disp(sprintf('g1=%d,field%d,chi2tot=%.3f',g1,ifield,chi2tot));
end
chi2dat(ig1).g1=g1;
chi2dat(ig1).chi2_arr=chi2_arr;
chi2dat(ig1).chi2tot=chi2tot;
end

%%
chi2dat_1st=chi2dat;
clear chi2dat
%% finner grid in G1
g1_arr=4:0.05:5.5;

for ig1=1:numel(g1_arr)
g1=-g1_arr(ig1);
chi2_arr=zeros(5,21);    
for ifield=4:8
mask=halfdat(ifield).bigmask;mkk=halfdat(ifield).wMkk;
rawdmap=halfdat(ifield).rawd;
nsimps=noise_sim_ps(flight,inst,ifield,mask,mkk,1,g1);
[psflight,~,cl2dflight]=get_Cl(rawdmap,mask,mkk,pixscale,nsimps.weight);

nmean=mean(nsimps.nCl_arr.*cal.*cal);nstd=std(nsimps.nCl_arr.*cal.*cal);
chi2=((psflight-nmean).^2)./(nstd.^2);
chi2_arr(ifield-3,:)=chi2;
chi2tot=sum(chi2(end-4:end));
disp(sprintf('g1=%d,field%d,chi2tot=%.3f',g1,ifield,chi2tot));
end
chi2dat(ig1).g1=g1;
chi2dat(ig1).chi2_arr=chi2_arr;
chi2dat(ig1).chi2tot=chi2tot;
end
chi2dat_2nd=chi2dat;
%%
save(strcat(savedir,'TM',num2str(inst),'_chi2data_1st'),'chi2dat_1st');
save(strcat(savedir,'TM',num2str(inst),'_chi2data_2nd'),'chi2dat_2nd');
%% plot chi2 
g1plot=zeros(1,50);
field4y=zeros(1,50);
field5y=zeros(1,50);
field6y=zeros(1,50);
field7y=zeros(1,50);
field8y=zeros(1,50);

for i=1:numel(chi2dat_1st)
g1plot(i)=chi2dat_1st(i).g1;
c=chi2dat_1st(i).chi2_arr;
field4y(i)=sum(c(1,end-4:end));
field5y(i)=sum(c(2,end-4:end));
field6y(i)=sum(c(3,end-4:end));
field7y(i)=sum(c(4,end-4:end));
field8y(i)=sum(c(5,end-4:end));
end

for i=1:numel(chi2dat_2nd)-2
g1plot(i+numel(chi2dat_1st))=chi2dat_2nd(i+1).g1;
c=chi2dat_2nd(i+1).chi2_arr;
field4y(i+numel(chi2dat_1st))=sum(c(1,end-4:end));
field5y(i+numel(chi2dat_1st))=sum(c(2,end-4:end));
field6y(i+numel(chi2dat_1st))=sum(c(3,end-4:end));
field7y(i+numel(chi2dat_1st))=sum(c(4,end-4:end));
field8y(i+numel(chi2dat_1st))=sum(c(5,end-4:end));
end

[g1plot,I]=sort(g1plot);
field4y=field4y(I);
field5y=field5y(I);
field6y=field6y(I);
field7y=field7y(I);
field8y=field8y(I);

%g1 4.5 is double counted
%g1plot(15)=[];field4y(15)=[];field5y(15)=[];field6y(15)=[];field7y(15)=[];
%field8y(15)=[];
bad=[15,16,20,22,24,26,28,30,32,34,36,37,39];
g1plot(bad)=[];field4y(bad)=[];field5y(bad)=[];
field6y(bad)=[];field7y(bad)=[];field8y(bad)=[];

fieldty=field4y+field5y+field6y+field7y+field8y;
%%
figure
plot(g1plot,fieldty,'-o');hold on
%plot(g1plot,field5y,'-o');
%plot(g1plot,field6y,'-o');
%plot(g1plot,field7y,'-o');
%plot(g1plot,field8y,'-o');
%plot(g1plot,fieldty,'-o');

minind=find(fieldty==min(fieldty));
mchi=fieldty(minind);pte=1-chi2cdf(mchi,24);
axis([-5.5 -3.5 0 1000])
xlabel('$-G1$','interpreter','latex','fontsize',18)
ylabel('$\chi^2$','interpreter','latex','fontsize',18)
title(sprintf('TM%d, %d, g1=%.2f, chi^2 = %.3f, PTE=%1.3f'...
            ,inst,flight,g1plot(minind),mchi,pte))
imname=sprintf('%sTM%dchi2tot',savedir,inst);
print(imname,'-dpng');%close
%% run again with fitted correct g1 diff
g1=-4.75;
chi2_arr=zeros(5,21);    
for ifield=4%4:8
mask=halfdat(ifield).bigmask;mkk=halfdat(ifield).wMkk;
rawdmap=halfdat(ifield).rawd;
nsimps=noise_sim_ps(flight,inst,ifield,mask,mkk,1,g1);
[psflight,~,cl2dflight,~,~,dpsflight]=get_Cl...
                (rawdmap,mask,mkk,pixscale,nsimps.weight);

nmean=mean(nsimps.nCl_arr.*cal.*cal);nstd=std(nsimps.nCl_arr.*cal.*cal);
chi2=((psflight-nmean).^2)./(nstd.^2);
chi2_arr(ifield-3,:)=chi2;
chi2tot=sum(chi2(end-4:end));
disp(sprintf('g1=%d,field%d,chi2tot=%.3f',g1,ifield,chi2tot));


l=nsimps.l;
rnCl_arr=nsimps.rnCl_arr.*cal.*cal;
phCl_arr=nsimps.phCl_arr.*cal.*cal;
nCl_arr=nsimps.nCl_arr.*cal.*cal;

loglog(l,l.*(l+1).*psflight./2./pi,'k.');hold on
errorbar(l,l.*(l+1).*psflight./2./pi,l.*(l+1).*dpsflight./2./pi,'.k',...
                                                    'markersize',20)
y1=(l.*(l+1).*(prctile(nCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(nCl_arr,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on
xlim([1e2,2e5]);ylim([1e-1,1e4]);
title(halfdat(ifield).name);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)

imname=sprintf('%sTM%d_i%d',savedir,inst,ifield);
print(imname,'-dpng');close
end

%% run again with fitted correct g1 full 

if inst==1;g1=-4.4;elseif inst==2;g1=-4.75;end
for ifield=4:8
mask=alldat(ifield).bigmask;mkk=alldat(ifield).wMkk;
calmap=alldat(ifield).calmap5;
nsimps=noise_sim_ps(flight,inst,ifield,mask,mkk,0,g1);
[psflight,~,cl2dflight,~,~,dpsflight]=get_Cl...
                (calmap,mask,mkk,pixscale,nsimps.weight);

nmean=mean(nsimps.nCl_arr.*cal.*cal);nstd=std(nsimps.nCl_arr.*cal.*cal);
chi2=((psflight-nmean).^2)./(nstd.^2);
chi2tot=sum(chi2(end-4:end));
disp(sprintf('g1=%d,field%d,chi2tot=%.3f',g1,ifield,chi2tot));


l=nsimps.l;
rnCl_arr=nsimps.rnCl_arr.*cal.*cal;
phCl_arr=nsimps.phCl_arr.*cal.*cal;
nCl_arr=nsimps.nCl_arr.*cal.*cal;

loglog(l,l.*(l+1).*psflight./2./pi,'k.');hold on
errorbar(l,l.*(l+1).*psflight./2./pi,l.*(l+1).*dpsflight./2./pi,'.k',...
                                                    'markersize',20)
y1=(l.*(l+1).*(prctile(nCl_arr,16))./2./pi);
y2=(l.*(l+1).*(prctile(nCl_arr,84))./2./pi);
pltn=fill([l,flip(l)],[abs(y1),abs(flip(y2))],...
    [0.2,0.6,0.2],'facealpha',0.5,'EdgeColor','none');hold on
xlim([1e2,2e5]);ylim([1e-2,5e3]);
title(halfdat(ifield).name);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)

imname=sprintf('%sTM%d_i%dfull',savedir,inst,ifield);
print(imname,'-dpng');close
end




