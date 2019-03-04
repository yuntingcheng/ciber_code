function [stampercb,stamperps,hitmap]=...
    stackihl_mock_star(flight,inst,ifield,Nsrc,m_min,m_max,interp,dx,...
    cbmap,psmap,mask_inst,strmask,verbose)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assign random star 
%
%Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - Nsrc: # of fake stars to be stacked
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - dx: stamp size is 2*dx+1, dx is 10 times finer CIBER pix
% - cbmap:
% - psmap:
% - mask_inst:
% - strmask:
% - strnum:
% - verbose:1 or 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);
catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PanSTARRS/');
%psfdir=(strcat(mypaths.ciberdir,'doc/20171130_psfstack/forward_model/'));
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

sp=find(x_arr>0.5 & x_arr<1024.5 & y_arr>0.5 & y_arr<1024.5);

x_arr = x_arr(sp);
y_arr = y_arr(sp);

%%% count the center pix map
xround_arr=round(x_arr);
yround_arr=round(y_arr);

centnum_map = zeros(1024);
for i=1:numel(sp)
    centnum_map(xround_arr(i),yround_arr(i))=...
        centnum_map(xround_arr(i),yround_arr(i))+1;
end

%%% assign the mock src position and mag %%%
sp = find(centnum_map == 0);
[subx_arr,suby_arr]=ind2sub([1024,1024],datasample(sp, Nsrc,'Replace',false));
subx_arr = subx_arr + rand(size(subx_arr))-0.5;
suby_arr = suby_arr + rand(size(suby_arr))-0.5;
xsmall_arr = subx_arr.*10 - 4.5;
ysmall_arr = suby_arr.*10 - 4.5;
subm_arr = rand(Nsrc,1)*(m_max-m_min) + m_min;

sr = ((7./3600.0)*(pi/180.0)).^2;

if interp==0
    lambdaeff=0.9633;
end

if interp==1
   if inst==1
       lambdaeff=1.05;
       ydat = squeeze(M(:,9)');
       Idat = squeeze(M(:,11)');
       sp = find(ydat>m_min & ydat<m_max);
       dm_dat = Idat(sp) - ydat(sp);
   else
       lambdaeff=1.79;
       ydat = squeeze(M(:,9)');
       Hdat = squeeze(M(:,12)');
       sp = find(ydat>m_min & ydat<m_max);
       dm_dat = Hdat(sp) - ydat(sp);
   end
   subm_arr = subm_arr + datasample(dm_dat,Nsrc)';   
end

subI_arr=3631*10.^(-subm_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);

%%% get psf params %%%
% psfparfile = strcat(psfdir,'TM',num2str(inst),'_',dt.name,'_bestparam.txt');
% params = csvread(psfparfile);
% A=params(1);
% B=params(2);
% sig=params(3);
% r0=params(4);
% alpha=params(5);

beta=fitpsfdat(ifield).psfmodel.beta_best;
rc=fitpsfdat(ifield).psfmodel.rc_best;
norm=fitpsfdat(ifield).psfmodel.norm;

stampercb=zeros(2*dx+1);
stamperps=zeros(2*dx+1);
hitmap=zeros(2*dx+1);

Nlarge = 10240+300+300;
radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
%Imap_large = A*exp(-radmap.^2./2./(sig)^2)+B./(1+(radmap./(r0)).^alpha);
Imap_large = norm .* (1 + (radmap/rc).^2).^(-3.*beta./2);

%%% sig clip params %%%
nbins = 25;
iter_clip = 3;
sig = 3;

print_count = 0;
for i=1:Nsrc    
    xi = round(xsmall_arr(i));
    yi = round(ysmall_arr(i));
    dx1 = Nlarge + 1 - xi;
    dy1 = Nlarge + 1 - yi;
    Imap = Imap_large(1+dx1:10240+dx1,1+dy1:10240+dy1).*subI_arr(i);
    srcmap = rebin_map_coarse(Imap,10).*100;

    cbmapi = (cbmap+srcmap).*strmask.*mask_inst;
    psmapi = (psmap+srcmap).*strmask.*mask_inst;
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');
    %%% rebin to 10x finer map
    stampcb0 = imresize(cbmapi,10,'nearest');
    stampps0 = imresize(psmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stampcb0 = stampcb0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stampps0 = stampps0(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    %%% sigma clip
    rmin = 20; % don't clip within rmin subpixels
    sig_clip_mask = stamp_clip(stampcb0,dx+1,dx+1,nbins,sig,iter_clip,rmin);
    stampcb = stampcb0 .* sig_clip_mask;
    stampps = stampps0 .* sig_clip_mask;
    maskstamp = sig_clip_mask;
    %%% stack
    stampercb=stampercb+stampcb;
    stamperps=stamperps+stampps;
    hitmap=hitmap+maskstamp;

    %%% print
    if verbose>0
        disp(sprintf('stack %d/%d src between %.1f<m<%.1f '...
            ,i,Nsrc,m_min,m_max));         
    end
 

%     setwinsize(gcf,1000,400)
%     subplot(121)
%     imageclip(stampercb);
%     subplot(122)
%     imageclip(stamperps);
%     title(i)
%     drawnow
end
    disp(sprintf('stack %d src between %.1f<m<%.1f '...
        ,Nsrc,m_min,m_max));

return
