flight=40030;
pixsizecb=7;
inst=2;
ifield=8;
savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/ukpurity/';

quad_arr=['A','B','C','D'];
dt=get_dark_times(flight,inst,ifield);

%SPHEREx
pixsize=6.2;
R=41.5;
%CDIM
%pixsize=1;
%=500;
%% get catalog to estimate number density
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/UKIDSS/');

iquad=1;
%%% read cat data %%%
quad=quad_arr(iquad);
catfile=strcat(catdir,dt.name,'_',quad,'_uk.txt');
M = csvread(catfile,1);
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end

if ifield==8
    cls_arr=squeeze(M(:,11)');
    ps_arr=squeeze(M(:,12)');
    pg_arr=squeeze(M(:,13)');
    pn_arr=squeeze(M(:,14)');
    psat_arr=squeeze(M(:,15)');
else
    cls_arr=squeeze(M(:,12)');
    ps_arr=squeeze(M(:,13)');
    pg_arr=squeeze(M(:,14)');
    pn_arr=squeeze(M(:,15)');
    psat_arr=squeeze(M(:,16)');    
end

figure
setwinsize(gcf,1000,400)
sp=find(pn_arr<0.04);
muse_arr=m_arr(sp);
binedge=5.5:24.5;
mbin_arr=(binedge(2:end)+binedge(1:end-1))/2;
Fbin_arr=3631.*10.^(mbin_arr/-2.5); %[Jy]
[Nbin_arr,edges] = histcounts(muse_arr,binedge);
Nbinuk_arr=Nbin_arr.*4./(1024)^2.*(pixsize/pixsizecb).^2;
subplot(1,2,1)
semilogy(mbin_arr,Nbinuk_arr);hold on
subplot(1,2,2)
semilogy(mbin_arr,Nbinuk_arr.*Fbin_arr.^2);hold on


catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PSC/');

iquad=1;
%%% read cat data %%%
quad=quad_arr(iquad);
catfile=strcat(catdir,dt.name,'_',quad,'_2m.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');
if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end


muse_arr=m_arr;
[Nbin_arr,edges] = histcounts(muse_arr,binedge);
Nbin2m_arr=Nbin_arr.*4./(1024)^2.*(pixsize/pixsizecb).^2;
subplot(1,2,1)
semilogy(mbin_arr,Nbin2m_arr);hold on
ylabel('$dN/dm$','interpreter','latex','fontsize',20)
xlabel('$m_{AB}$','interpreter','latex','fontsize',20);

subplot(1,2,2)
semilogy(mbin_arr,Nbin2m_arr.*Fbin_arr.^2);
ylabel('$S^2dN/dm [Jy^2]$','interpreter','latex','fontsize',20)
xlabel('$m_{AB}$','interpreter','latex','fontsize',20);


Nbin_arr=Nbin2m_arr;
Nbin_arr(find(mbin_arr>17))=Nbinuk_arr(find(mbin_arr>17));
subplot(1,2,1)
semilogy(mbin_arr,Nbin_arr,'o');

subplot(1,2,2)
semilogy(mbin_arr,Nbin_arr.*Fbin_arr.^2,'o');

%%
figure
setwinsize(gcf,1000,400)

Sint_arr=zeros(1,numel(mbin_arr));
Nint_arr=zeros(1,numel(mbin_arr));
for i=1:numel(mbin_arr)
    Sint_arr(i)=sum(Nbin_arr(i:end).*Fbin_arr(i:end).^2);
    Nint_arr(i)=sum(Nbin_arr(1:i));
end


subplot(1,2,1)
semilogy(mbin_arr,sqrt(Sint_arr),'-o');
ylabel('$\sqrt{\int dm S^2dN/dm (>m)} [Jy]$','interpreter','latex','fontsize',20)
xlabel('$m_{AB}$','interpreter','latex','fontsize',20);

subplot(1,2,2)
semilogy(mbin_arr,Nint_arr,'-o');
ylabel('$N(<m)$','interpreter','latex','fontsize',20)
xlabel('$m_{AB}$','interpreter','latex','fontsize',20);

%%
% H alpha
nu_e=3e8/0.656e-6;
z=2.23;
DL=1.24e4; %[Mpc/h]
nu_obs=nu_e*(1+z);
L_st=10^42.87;%[erg/s]
L_st=L_st*1e-7/3.828e26;%[M_sun]


sig2=sum(Nbin_arr.*Fbin_arr.^2);
sig=sqrt(sig2); %[Jy]
sig=sig*nu_obs/R*1e-26;%[W/m^2]
sig=sig*(3.1e22)^2/3.828e26; % [Msun/Mpc^2]
sigL=sig*4*pi*(DL/0.7)^2; %[Msun]

Lbin_arr=Fbin_arr.*nu_obs./R.*1e-26;
Lbin_arr=Lbin_arr.*(3.1e22).^2./3.828e26;
Lbin_arr=Lbin_arr.*4.*pi.*(DL/0.7).^2;

figure
semilogy(mbin_arr,Lbin_arr./L_st);hold on
semilogy(mbin_arr,ones(size(mbin_arr)).*sigL./L_st);
ylabel('$\sigma_L/L_*$','interpreter','latex','fontsize',20)
xlabel('$m_{AB}$','interpreter','latex','fontsize',20);
