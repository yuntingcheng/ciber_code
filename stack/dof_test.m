%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian, no masking
dx = 1200;
Nsim = 100;
Nstack = 100;
flight = 40030;
mypaths = get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/DoF_test/'));

for j=1:Nsim
map = randn(1024);
mask = ones(1024);
subx_arr=(dx/10 + 1) + (1023 - dx/5) * rand(1,Nstack);
suby_arr=(dx/10 + 1) + (1023 - dx/5) * rand(1,Nstack);
stamper=zeros(2*dx+1);
stamper2=zeros(2*dx+1);
maskstamper=zeros(2*dx+1);

for i=1:Nstack
    disp(sprintf('sim %d, stack %d', j, i));
    mapi = map.*mask;
    %%% rebin to 10x finer map
    stamp = imresize(mapi,10,'nearest');
    stamp2 = randn(10240);
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5);
    ycent = round(suby_arr(i)*10-4.5);
    stamp = stamp(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stamp2 = stamp2(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp=zeros(size(stamp));
    maskstamp(find(stamp))=1;
    %%% stack
    if mod(i,4)==0
        stamper=stamper+imrotate(stamp, 0);
        stamper2=stamper2+imrotate(stamp2, 0);
        maskstamper=maskstamper+imrotate(maskstamp, 0);
    elseif mod(i,4)==1
        stamper=stamper+imrotate(stamp, 90);
        stamper2=stamper2+imrotate(stamp2, 90);
        maskstamper=maskstamper+imrotate(maskstamp, 90);
    elseif mod(i,4)==2
        stamper=stamper+imrotate(stamp, 180);
        stamper2=stamper2+imrotate(stamp2, 180);
        maskstamper=maskstamper+imrotate(maskstamp, 180);
    elseif mod(i,4)==3
        stamper=stamper+imrotate(stamp, 270);
        stamper2=stamper2+imrotate(stamp2, 270);
        maskstamper=maskstamper+imrotate(maskstamp, 270);
    end 
end
profile = radial_prof(stamper./maskstamper,ones(2*dx+1),...
    dx+1,dx+1,1,25);    
Mprof1(j,:) = profile.prof;
Merr1(j,:) = profile.err;
profile = radial_prof(stamper2./maskstamper,ones(2*dx+1),...
    dx+1,dx+1,1,25);
Mprof2(j,:) = profile.prof;
Merr2(j,:) = profile.err;
end
r_arr = profile.r*0.7;
%%
figure
setwinsize(gcf,800,300)
subplot(1,2,1)
loglog(r_arr, std(Mprof2),'b','linewidth', 3);hold on
loglog(r_arr, Merr2(1,:),'r');
legend({'\sigma_{sims}','\sigma_{annuli} / sqrt(N)'},...
    'Location','southwest')

for i=2:Nsim
    loglog(r_arr, Merr2(i,:),'r');
end
loglog(r_arr, std(Mprof2),'b','linewidth', 3);
title('no subpixel stacking','fontsize',20);
xlim([4e-1,1e3])
ylim([5e-5,1e-1])
xlabel('arcsec', 'fontsize',20);
ylabel('\sigma_{stack}','fontsize',20);

subplot(1,2,2)
loglog(r_arr, std(Mprof1),'b','linewidth', 3);hold on
loglog(r_arr, Merr1(1,:),'r');
loglog(r_arr, std(Mprof1)./sqrt(100),'b--', 'linewidth', 3);hold on
legend({'\sigma_{sims}','\sigma_{annuli} / sqrt(N)','\sigma_{sims}/sqrt(100)'},...
    'Location','southwest')
for i=2:Nsim
    loglog(r_arr, Merr1(i,:),'r');
end
loglog(r_arr, std(Mprof1),'b','linewidth', 3);
loglog(r_arr, std(Mprof1)./sqrt(100),'b--', 'linewidth', 3);
title('10x10 subpixel stacking','fontsize',20);
xlim([4e-1,1e3])
ylim([5e-5,1e-1])
xlabel('arcsec', 'fontsize',20);
ylabel('<I_{stack}>','fontsize',20);

savename=strcat(pltsavedir,'Gaussian_no_mask');
print(savename,'-dpng');%close
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gaussian, with mask
dx = 1200;
Nsim = 100;
Nstack = 100;
flight = 40030;
inst = 1;
ifield = 8;
mypaths = get_paths(flight);
pltsavedir=(strcat(mypaths.ciberdir,'doc/20171018_stackihl/DoF_test/'));

dt=get_dark_times(flight,inst,ifield);
loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(sprintf('%sstackmapdat',loaddir),'stackmapdat');
mask = stackmapdat(ifield).mask_inst_clip.*stackmapdat(ifield).strmask;
clear stackmapdat
%%
for j=1:Nsim
map = randn(1024);
subx_arr=(dx/10 + 1) + (1023 - dx/5) * rand(1,Nstack);
suby_arr=(dx/10 + 1) + (1023 - dx/5) * rand(1,Nstack);
stamper=zeros(2*dx+1);
maskstamper=zeros(2*dx+1);

for i=1:Nstack
    disp(sprintf('sim %d, stack %d', j, i));
    mapi = map.*mask;
    %%% rebin to 10x finer map
    stamp = imresize(mapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5);
    ycent = round(suby_arr(i)*10-4.5);
    stamp = stamp(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp=zeros(size(stamp));
    maskstamp(find(stamp))=1;
    %%% stack
    if mod(i,4)==0
        stamper=stamper+imrotate(stamp, 0);
        maskstamper=maskstamper+imrotate(maskstamp, 0);
    elseif mod(i,4)==1
        stamper=stamper+imrotate(stamp, 90);
        maskstamper=maskstamper+imrotate(maskstamp, 90);
    elseif mod(i,4)==2
        stamper=stamper+imrotate(stamp, 180);
        maskstamper=maskstamper+imrotate(maskstamp, 180);
    elseif mod(i,4)==3
        stamper=stamper+imrotate(stamp, 270);
        maskstamper=maskstamper+imrotate(maskstamp, 270);
    end 
end
profile = radial_prof(stamper./maskstamper,ones(2*dx+1),...
    dx+1,dx+1,1,25,'weight',maskstamper);    
Mprof(j,:) = profile.prof;
Merr(j,:) = profile.err;
end
r_arr = profile.r*0.7;
%%
figure
setwinsize(gcf,800,300)
subplot(1,2,1)
imageclip(mask.*map);

subplot(1,2,2)
loglog(r_arr, std(Mprof),'b','linewidth', 3);hold on
loglog(r_arr, Merr(1,:),'r');
loglog(r_arr, std(Mprof)./sqrt(100),'b--', 'linewidth', 3);
legend({'\sigma_{sims}','\sigma_{annuli} / sqrt(N)','\sigma_{sims}/sqrt(100)'},...
    'Location','southwest')

for i=2:Nsim
    loglog(r_arr, Merr(i,:),'r');
end
loglog(r_arr, std(Mprof),'b','linewidth', 3);
loglog(r_arr, std(Mprof)./sqrt(100),'b--', 'linewidth', 3);
xlim([4e-1,1e3])
ylim([1e-4,1e0])
xlabel('arcsec', 'fontsize',20);
ylabel('\sigma_{stack}','fontsize',20);

h = suptitle('masked Gaussian random field');
set(h,'FontSize',20);

savename=strcat(pltsavedir,'Gaussian_mask');
%print(savename,'-dpng');%close