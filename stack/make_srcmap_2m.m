function map = make_srcmap_2m(flight,inst,ifield,m_min,m_max,interp,PSmatch)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the sim src map from 2MASS
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8
% - m_min: min masking magnitude
% - m_max: max masking magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/PSC/');
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.csv');

M = csvread(catfile,1);

if inst==1
    x_arr=squeeze(M(:,4)');
    y_arr=squeeze(M(:,3)');
else
    x_arr=squeeze(M(:,6)');
    y_arr=squeeze(M(:,5)');    
end
x_arr=x_arr+1;
y_arr=y_arr+1;

% y band
m_arr=squeeze(M(:,12)');
sr = ((7./3600.0)*(pi/180.0)).^2;

% use y band mag
if interp == 0
    lambdaeff=0.9633;
    I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
end

% use linear interpolated magnitude
if interp == 1
    if inst==1
        mcb_arr = squeeze(M(:,10)');
        lambdaeff=1.05;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    else
        mcb_arr = squeeze(M(:,11)');
        lambdaeff=1.79;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    end
end

if interp == 2
    if inst==1
        mcb_arr = squeeze(M(:,10)');
        lambdaeff=1.05;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    else
        mcb_arr = squeeze(M(:,11)');
        lambdaeff=1.79;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    end
    m_arr = mcb_arr;
end

%%% select cat data %%%
sp=find(m_arr<=m_max & m_arr>m_min ...
     & x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);

if PSmatch~=0
    match_arr = squeeze(M(:,11)');
    sp=find(m_arr<=m_max & m_arr>m_min & match_arr==0 ...
     & x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);
end

subI_arr=I_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

%%% get x,y coord in small grid %%%
xsmall_arr = subx_arr.*10 - 4.5;
ysmall_arr = suby_arr.*10 - 4.5;

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

%%% make srcmap %%%
Nsrc = numel(subI_arr);
srcmap=zeros(10240);
disp(sprintf('make srcmap TM%d %s, mrange=(%d,%d), %d srcs',...
    inst,dt.name,m_min,m_max,Nsrc));
if Nsrc>0
    Nlarge = 10240+300+300;
    radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
    %Imap_large = A*exp(-radmap.^2./2./(sig)^2)+B./(1+(radmap./(r0)).^alpha);
    Imap_large = norm .* (1 + (radmap./rc).^2).^(-3.*beta./2);
    if Nsrc>20
        idx_print=floor(Nsrc/20) * (1:20);
    else
        idx_print=[];
    end

    print_count=0;  
    for i=1:Nsrc
        xi = round(xsmall_arr(i));
        yi = round(ysmall_arr(i));
        dx = Nlarge + 1 - xi;
        dy = Nlarge + 1 - yi;
        Imap = Imap_large(1+dx:10240+dx,1+dy:10240+dy).*subI_arr(i);
        srcmap = srcmap + Imap;

        if ismember(i,idx_print)
            print_count = print_count + 1;
            disp(sprintf('make srcmap %d %% sources',print_count*5));
        end
    end
end
map = rebin_map_coarse(srcmap,10);
map = map.*100;
end