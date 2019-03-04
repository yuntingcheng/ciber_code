flight=40030;
inst=1;
mypaths=get_paths(flight);

psfdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf/yt/inst',...
        num2str(inst),'/j0_14/');
% old I band psf from PMK
%psfdir='/Users/ytcheng/ciber/data/iband_psf_quads/';

savedir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/blfit/TM',...
    num2str(inst),'/');
load(strcat(savedir,'Cldat'),'Cldat');
%% plot bl
pixsize=0.7;
for ifield=4:8

    bl_arr=zeros(4,21);
    dt=get_dark_times(flight,inst,ifield);
    name=dt.name;
    
    a=fitsread(strcat(psfdir,name,'_A.fits'));
    b=fitsread(strcat(psfdir,name,'_B.fits'));
    c=fitsread(strcat(psfdir,name,'_C.fits'));
    d=fitsread(strcat(psfdir,name,'_D.fits'));
    
    ah=fitsread(strcat(psfdir,name,'_Ahitmap.fits'));
    bh=fitsread(strcat(psfdir,name,'_Bhitmap.fits'));
    ch=fitsread(strcat(psfdir,name,'_Chitmap.fits'));
    dh=fitsread(strcat(psfdir,name,'_Dhitmap.fits'));

    psfmap=(a.*ah+b.*bh+c.*ch+d.*dh)./(ah+bh+ch+dh);

    [bla,lfine]=get_angular_spec(a,a,pixsize,'nbins',30);
    [blb]=get_angular_spec(b,b,pixsize,'nbins',30);
    [blc]=get_angular_spec(c,c,pixsize,'nbins',30);
    [bld]=get_angular_spec(d,d,pixsize,'nbins',30);
    [blt]=get_angular_spec(psfmap,psfmap,pixsize,'nbins',30);
    
    lfine=lfine(find(bla));
    bla=bla(find(bla));blb=blb(find(blb));blc=blc(find(blc));bld=bld(find(bld));
    blt=blt(find(blt));
    
    figure
    setwinsize(gcf,1000,600)
    subplot(2,3,1)
    loglog(lfine,bla,'mo-');hold on
    loglog(lfine,blb,'ro-');hold on
    loglog(lfine,blc,'go-');hold on
    loglog(lfine,bld,'bo-');hold on
    loglog(lfine,blt,'k+-');hold on
    xlim([1e3,2e6])
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_\ell$','interpreter','latex','fontsize',18)
    title('PSF map power spectrum')
    
    subplot(2,3,4)
    semilogx(lfine,bla,'mo-');hold on
    loglog(lfine,blb,'ro-');hold on
    loglog(lfine,blc,'go-');hold on
    loglog(lfine,bld,'bo-');hold on
    loglog(lfine,blt,'k+-');hold on
    xlim([1e3,2e6])
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_\ell$','interpreter','latex','fontsize',18)

    bla(1:2)=mean(bla(3:4));blb(1:2)=mean(blb(3:4));
    blc(1:2)=mean(blc(3:4));bld(1:2)=mean(bld(3:4)); 
    blt(1:2)=mean(blt(3:4));
    
    blaf=bla./mean(bla(3:4));blbf=blb./mean(blb(3:4));
    blcf=blc./mean(blc(3:4));bldf=bld./mean(bld(3:4));
    bltf=blt./mean(blt(3:4));
    
    subplot(2,3,2)
    loglog(lfine,blaf,'mo-');hold on
    loglog(lfine,blbf,'ro-');hold on
    loglog(lfine,blcf,'go-');hold on
    loglog(lfine,bldf,'bo-');hold on
    loglog(lfine,bltf,'k+-');hold on
    xlim([1e3,2e6])
    ylim([1e-7,1.5])
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_\ell$','interpreter','latex','fontsize',18)
    title('normalized bl')
        
    subplot(2,3,5)
    semilogx(lfine,blaf,'mo-');hold on
    loglog(lfine,blbf,'ro-');hold on
    loglog(lfine,blcf,'go-');hold on
    loglog(lfine,bldf,'bo-');hold on
    loglog(lfine,bltf,'k+-');hold on
    xlim([1e3,2e6])
    ylim([0,1.1])
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_\ell$','interpreter','latex','fontsize',18)
    
    l=Cldat(ifield).l;
    l1=[l(1) lfine];
    
    bl1=[1 blaf];
    logbl=interp1(log10(l1),log10(bl1),log10(l),'pchip','extrap');
    bla=10.^logbl;
    bl_arr(1,:)=bla;
    
    bl1=[1 blbf];
    logbl=interp1(log10(l1),log10(bl1),log10(l),'pchip','extrap');
    blb=10.^logbl;
    bl_arr(2,:)=blb;
    
    bl1=[1 blcf];
    logbl=interp1(log10(l1),log10(bl1),log10(l),'pchip','extrap');
    blc=10.^logbl;
    bl_arr(3,:)=blc;
    
    bl1=[1 bldf];
    logbl=interp1(log10(l1),log10(bl1),log10(l),'pchip','extrap');
    bld=10.^logbl;
    bl_arr(4,:)=bld;

    bl1=[1 bltf];
    logbl=interp1(log10(l1),log10(bl1),log10(l),'pchip','extrap');
    blt=10.^logbl;

    subplot(2,3,3)
    loglog(l,bla,'m-','DisplayName','quad 1');hold on
    loglog(l,blb,'r-','DisplayName','quad 2');hold on
    loglog(l,blc,'g-','DisplayName','quad 3');hold on
    loglog(l,bld,'b-','DisplayName','quad 4');hold on
    loglog(l,blt,'k+-','DisplayName','stack quad');hold on
    loglog(l,Cldat(ifield).blmz,'k.-','DisplayName','bl MZ');hold on
    xlim([1e2,2e5])
    ylim([1e-3,1.5])
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_\ell$','interpreter','latex','fontsize',18)
    h=legend('show','Location','southwest');
    set(h,'fontsize',10)
    legend boxoff
    title('bl extrapolation')
    
    subplot(2,3,6)
    semilogx(l,bla,'m-');hold on
    loglog(l,blb,'r-');hold on
    loglog(l,blc,'g-');hold on
    loglog(l,bld,'b-');hold on
    loglog(l,blt,'k+-');hold on
    loglog(l,Cldat(ifield).blmz,'k.-');hold on
    
    xlim([1e2,2e5])
    ylim([0,1.1])
    xlabel('$\ell$','interpreter','latex','fontsize',18)
    ylabel('$b_\ell$','interpreter','latex','fontsize',18)

    drawnow
    savename=strcat(savedir,'bl_i',num2str(ifield));
    print(savename,'-dpng');close    

    Cldat(ifield).blmax=max(bl_arr);
    Cldat(ifield).blmin=min(bl_arr);
    Cldat(ifield).blt=blt;
    
end
save(strcat(savedir,'Cldat'),'Cldat');
%% plot flight PS no scale RN or ph
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
G1=-2.9;

for ifield=4:8
fig=figure;
dt=get_dark_times(flight,inst,ifield);
l=Cldat(ifield).l;
fCl=Cldat(ifield).flightCl.*cal.^2;
rnClFF=Cldat(ifield).rnClFF.*cal.^2;
phClFF=Cldat(ifield).phClFF./(-G1).*cal.^2;
drnClFF=Cldat(ifield).drnClFF.*cal.^2;
dphClFF=Cldat(ifield).dphClFF./(-G1).*cal.^2;
rnCl=Cldat(ifield).rnCl.*cal.^2;
phCl=Cldat(ifield).phCl./(-G1).*cal.^2;
drnCl=Cldat(ifield).drnCl.*cal.^2;
dphCl=Cldat(ifield).dphCl./(-G1).*cal.^2;

Cldb=fCl-rnCl-phCl-rnClFF-phClFF;
dCldb=sqrt(drnCl.^2+dphCl.^2+drnClFF.^2+dphClFF.^2);

%%%%% (Raw flight -FFbias - noise)/bl mean quad%%%%%
yvalue=(l.*(l+1).*Cldb./Cldat(ifield).blt./2./pi);
ebar=(l.*(l+1).*dCldb./Cldat(ifield).blt./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii) & yvalue(ii)>0
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    if yvalue(ii)<0 & yvalue(ii)+ebar(ii)>0
        upper=yvalue(ii)+ebar(ii);
        yvalue(ii)=1e-8;
        ebar(ii)=upper-yvalue(ii);
    end
end
errorbar(l,yvalue,ebarlow,ebar,'ko',...
'DisplayName','stack quad','markersize',3);hold on

%%%%% (Raw flight -FFbias - noise)/blmz%%%%%
yvalue=(l.*(l+1).*Cldb./Cldat(ifield).blmz./2./pi);
ebar=(l.*(l+1).*dCldb./Cldat(ifield).blmz./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii) & yvalue(ii)>0
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    if yvalue(ii)<0 & yvalue(ii)+ebar(ii)>0
        upper=yvalue(ii)+ebar(ii);
        yvalue(ii)=1e-8;
        ebar(ii)=upper-yvalue(ii);
    end
end
errorbar(l.*1.05,yvalue,ebarlow,ebar,'mo',...
'DisplayName','bl MZ','markersize',3);hold on
%%%%% (Raw flight -FFbias - noise)/blmz%%%%%
yvalue=(l.*(l+1).*Cldb./Cldat(ifield).blmax./2./pi);
ebar=(l.*(l+1).*dCldb./Cldat(ifield).blmax./2./pi);
ebarlow=ebar;
for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii) & yvalue(ii)>0
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    if yvalue(ii)<0 & yvalue(ii)+ebar(ii)>0
        upper=yvalue(ii)+ebar(ii);
        yvalue(ii)=1e-8;
        ebar(ii)=upper-yvalue(ii);
    end
end
errorbar(l.*1.1,yvalue,ebarlow,ebar,'bo',...
'DisplayName','bl max quad','markersize',3);hold on

%%%%% (Raw flight -FFbias - noise)/blmin%%%%%
yvalue=(l.*(l+1).*Cldb./Cldat(ifield).blmin./2./pi);
ebar=(l.*(l+1).*dCldb./Cldat(ifield).blmin./2./pi);
ebarlow=ebar;

for ii=1:numel(l)
    if yvalue(ii)<=ebar(ii) & yvalue(ii)>0
        ebarlow(ii)=yvalue(ii)-1e-8;
    end
    if yvalue(ii)<0 & yvalue(ii)+ebar(ii)>0
        upper=yvalue(ii)+ebar(ii);
        yvalue(ii)=1e-8;
        ebar(ii)=upper-yvalue(ii);
    end
end
errorbar(l.*1.15,yvalue,ebarlow,ebar,'ro',...
'DisplayName','bl min quad','markersize',3);hold on


xlim([1e2,2e5]);ylim([1e0,1e4]);
ax = get(fig,'CurrentAxes');
set(ax,'XScale','log','YScale','log');

title(strcat(dt.name));
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff

savename=strcat(savedir,'unscale_i',num2str(ifield));
print(savename,'-dpng');close    

end
%% fix RN, change G1
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
for ifield=4:8
Arn=1;
dt=get_dark_times(flight,inst,ifield);

figure
for G1=-3.9:0.2:-2.9
l=Cldat(ifield).l;
fCl=Cldat(ifield).flightCl;
rnClFF=Cldat(ifield).rnClFF.*Arn;
phClFF=Cldat(ifield).phClFF./(-G1);
drnClFF=Cldat(ifield).drnClFF.*Arn;
dphClFF=Cldat(ifield).dphClFF./(-G1);
rnCl=Cldat(ifield).rnCl.*Arn;
phCl=Cldat(ifield).phCl./(-G1);
drnCl=Cldat(ifield).drnCl.*Arn;
dphCl=Cldat(ifield).dphCl./(-G1);

bl=Cldat(ifield).blt;

Cldb=fCl-rnCl-phCl-rnClFF-phClFF;
dCldb=sqrt(drnCl.^2+dphCl.^2+drnClFF.^2+dphClFF.^2);

loglog(l,l.*(l+1).*Cldb./bl.*cal.^2./2./pi,'.-',...
    'DisplayName',strcat('G1=',num2str(G1)));hold on

end
norm=l.*(l+1).*Cldb./bl.*cal.^2./2./pi;norm=norm(16);
refline=l.*(l+1);refline=refline.*norm/refline(16);
loglog(l,refline,'k--','DisplayName','white slope');hold on
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)

xlim([1e2,2e5]);ylim([1e0,1e5]);
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
title(dt.name);

savename=strcat(savedir,'scaleG1_i',num2str(ifield));
print(savename,'-dpng');close    

end
%%
G1_arr=-3.9:0.2:-2.9;

for ifield=4:8
figure
setwinsize(gcf,1000,600)
count=0;
for Arn=0.8:0.2:1.4
count=count+1;
subplot(2,2,count)
for iG1=1:numel(G1_arr)
G1=G1_arr(iG1);
l=Cldat(ifield).l;
fCl=Cldat(ifield).flightCl;
rnClFF=Cldat(ifield).rnClFF.*Arn;
phClFF=Cldat(ifield).phClFF./(-G1);
drnClFF=Cldat(ifield).drnClFF.*Arn;
dphClFF=Cldat(ifield).dphClFF./(-G1);
rnCl=Cldat(ifield).rnCl.*Arn;
phCl=Cldat(ifield).phCl./(-G1);
drnCl=Cldat(ifield).drnCl.*Arn;
dphCl=Cldat(ifield).dphCl./(-G1);

bl=Cldat(ifield).blt;

Cldb=fCl-rnCl-phCl-rnClFF-phClFF;
dCldb=sqrt(drnCl.^2+dphCl.^2+drnClFF.^2+dphClFF.^2);


if iG1==2
e=errorbar(l(10:end).*0.95,Cldb(10:end)./bl(10:end),dCldb(10:end)./bl(10:end),...
    'DisplayName',strcat('G1=',num2str(G1)));hold on
e.Marker='.';e.Color=get_color(iG1);
elseif iG1==3
e=errorbar(l(10:end),Cldb(10:end)./bl(10:end),dCldb(10:end)./bl(10:end),...
    'DisplayName',strcat('G1=',num2str(G1)));hold on
e.Marker='.';e.Color=get_color(iG1);
elseif iG1==4
e=errorbar(l(10:end).*1.05,Cldb(10:end)./bl(10:end),dCldb(10:end)./bl(10:end),...
    'DisplayName',strcat('G1=',num2str(G1)));hold on
e.Marker='.';e.Color=get_color(iG1);
    
else
 semilogx(l(10:end),Cldb(10:end)./bl(10:end),'.-','color',get_color(iG1),...
    'DisplayName',strcat('G1=',num2str(G1)));hold on
end 
drawnow
end
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$C_\ell$','interpreter','latex','fontsize',18)

ylim([0,3e-10])
xlim([2e3,2e5])
h=legend('show','Location','northwest');
set(h,'fontsize',10)
legend boxoff
title(sprintf('Arn=%.2f',Arn));
end
savename=strcat(savedir,'scaleG1rn_i',num2str(ifield));
print(savename,'-dpng');close    
end