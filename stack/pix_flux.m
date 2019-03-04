%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%comparing CIBER and UKIDSS sim map flux pixel by pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
field=4;
quad_arr=['A','B','C','D'];
%%
masktot=zeros(1024);
cbmaptot=zeros(1024);
ukmaptot=zeros(1024);

figure
setwinsize(gcf,600,600)
for iquad=1:numel(quad_arr)
quad=quad_arr(iquad);

%%%%% load the mask and maps %%%%%
dt=get_dark_times(flight,inst,field);
outdir='/Users/ytcheng/ciber/doc/20170617_Stacking/maps/quadmasks/';
file=strcat(outdir,dt.name,'_inst',num2str(inst),'_mask_',quad,'.fits');
mask=fits_read(file);
outdir='/Users/ytcheng/ciber/doc/20170617_Stacking/maps/quadmaps/';
file=strcat(outdir,dt.name,'_inst',num2str(inst),'_map_',quad,'.fits');
cbmap=fits_read(file);
ukdir='/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/';
file=strcat(ukdir,'TM',num2str(inst),'/',dt.name,'_',quad,'_srcmapt.fits');
ukmap=fits_read(file);

%%%%% stick the whole map %%%%%
if iquad==1
    masktot(1:512,1:512)=mask;
    cbmaptot(1:512,1:512)=cbmap;
    ukmaptot(1:512,1:512)=ukmap;
elseif iquad==2
    masktot(513:1024,1:512)=mask;
    cbmaptot(513:1024,1:512)=cbmap;
    ukmaptot(513:1024,1:512)=ukmap;
elseif iquad==3
    masktot(1:512,513:1024)=mask;
    cbmaptot(1:512,513:1024)=cbmap;
    ukmaptot(1:512,513:1024)=ukmap;
else
    masktot(513:1024,513:1024)=mask;
    cbmaptot(513:1024,513:1024)=cbmap;
    ukmaptot(513:1024,513:1024)=ukmap;  
end
    
%%%%% start binning %%%%%
Iedge_arr=logspace(log10(30),log10(1e6),21);
Iedge_arr=[Iedge_arr,1e20];
x_arr=zeros(1,numel(Iedge_arr)-1);
ex_arr=zeros(1,numel(Iedge_arr)-1);
y_arr=zeros(1,numel(Iedge_arr)-1);
ey_arr=zeros(1,numel(Iedge_arr)-1);

for iI=1:numel(Iedge_arr)-1
    sp=find(ukmap>Iedge_arr(iI) & ukmap<=Iedge_arr(iI+1) & mask==1);
    cbuse=cbmap(sp);
    %{
    medcb=median(cbuse);
    sigcb=std(cbuse);
    cbgood=find(cbuse<medcb+3*sigcb);
    cbuse=cbuse(cbgood);
    %}
    x_arr(iI)=median(ukmap(sp));
    ex_arr(iI)=std(ukmap(sp))./sqrt(numel(sp));
    y_arr(iI)=median(cbuse);
    ey_arr(iI)=std(cbuse)./sqrt(numel(sp));
end

subplot(2,2,iquad)
errorbar(x_arr,y_arr,ey_arr,ey_arr,ex_arr,ex_arr,'k.');hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([10,1e7])
plot([10,1e7],[10,1e7]);
title(strcat(dt.name,'\_',quad),'fontsize',12)
xlabel('UKIDSS sim ($nW/m^2/sr$)','interpreter','latex','fontsize',18)
ylabel('CIBER data ($nW/m^2/sr$)','interpreter','latex','fontsize',18)
end
%%
figure
%setwinsize(gcf,600,600)
for iquad=1:numel(quad_arr)
quad=quad_arr(iquad);

%%%%% load the mask and maps %%%%%
dt=get_dark_times(flight,inst,field);
outdir='/Users/ytcheng/ciber/doc/20170617_Stacking/maps/quadmasks/';
file=strcat(outdir,dt.name,'_inst',num2str(inst),'_mask_',quad,'.fits');
mask=fits_read(file);
outdir='/Users/ytcheng/ciber/doc/20170617_Stacking/maps/quadmaps/';
file=strcat(outdir,dt.name,'_inst',num2str(inst),'_map_',quad,'.fits');
cbmap=fits_read(file);
ukdir='/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/';
file=strcat(ukdir,'TM',num2str(inst),'/',dt.name,'_',quad,'_srcmapt.fits');
ukmap=fits_read(file);
%%%%% start binning %%%%%
Iedge_arr=logspace(log10(30),log10(1e6),21);
Iedge_arr=[Iedge_arr,1e20];
x_arr=zeros(1,numel(Iedge_arr)-1);
ex_arr=zeros(1,numel(Iedge_arr)-1);
y_arr=zeros(1,numel(Iedge_arr)-1);
ey_arr=zeros(1,numel(Iedge_arr)-1);

for iI=1:numel(Iedge_arr)-1
    sp=find(cbmap>Iedge_arr(iI) & cbmap<=Iedge_arr(iI+1) & mask==1);
    ukuse=ukmap(sp);
    x_arr(iI)=median(cbmap(sp));
    ex_arr(iI)=std(cbmap(sp))./sqrt(numel(sp));
    y_arr(iI)=median(ukuse);
    ey_arr(iI)=std(ukuse)./sqrt(numel(sp));
end

%subplot(2,2,iquad)
errorbar(x_arr,y_arr,ey_arr,ey_arr,ex_arr,ex_arr,'.');hold on
set(gca,'xscale','log')
set(gca,'yscale','log')
xlim([10,1e7])
plot([10,1e7],[10,1e7]);
title(strcat(dt.name,'\_',quad),'fontsize',12)
xlabel('CIBER data ($nW/m^2/sr$)','interpreter','latex','fontsize',18)
ylabel('UKIDSS sim ($nW/m^2/sr$)','interpreter','latex','fontsize',18)
end
