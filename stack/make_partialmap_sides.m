function [srcmap_out, ihlmap_out,mask] = make_partialmap_sides...
    (flight,inst,ifield,m_min,m_max,interp,params,radius,frac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the sim src map PanSTARR
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - type: 1:gal, -1:star, 0:all, 2:undefined
% - m_min: min masking magnitude (PS y band)
% - m_max: max masking magnitude (PS y band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,'sides.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

m_arr=squeeze(M(:,4)');

sr = ((7./3600.0)*(pi/180.0)).^2;

% use y band mag
if interp == 0
    lambdaeff=0.9633;
    I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
end

% use linear interpolated magnitude
if interp == 1
    if inst==1
        mcb_arr = squeeze(M(:,5)');
        lambdaeff=1.05;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    else
        mcb_arr = squeeze(M(:,6)');
        lambdaeff=1.79;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    end
    
end

%%% select cat data %%%
sp_all=find(m_arr<=m_max & m_arr>m_min);

Nstack = round(frac * numel(sp_all));
randidx = randperm(numel(sp_all),Nstack);
sp = sp_all(randidx);

subI_arr=I_arr(sp);
subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

%%% get x,y coord in small grid %%%
xsmall_arr = subx_arr.*10 - 4.5;
ysmall_arr = suby_arr.*10 - 4.5;

%%%%%%%%%%%%%%%%%%%%%%%%%% make srcmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get psf params %%%
psfdir=(strcat(mypaths.ciberdir, 'doc/20171130_psfstack/forward_model/'));
psfparfile = strcat(psfdir,'TM',num2str(inst),'_',dt.name,'_bestparam.txt');
psfparams = csvread(psfparfile);
A=psfparams(1);
B=psfparams(2);
sig=psfparams(3);
r0=psfparams(4);
alpha=psfparams(5);

%%% make srcmap %%%
Nsrc = numel(subI_arr);
srcmap=zeros(7200);
disp(sprintf('make srcmap TM%d %s, mrange=(%d,%d), %d srcs',...
    inst,dt.name,m_min,m_max,Nsrc));
if Nsrc>0
    Nlarge =10240+300+300;
    radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
    Imap_large = A*exp(-radmap.^2./2./(sig)^2)+B./(1+(radmap./(r0)).^alpha);

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
        Imap = Imap_large(1+dx:7200+dx,1+dy:7200+dy).*subI_arr(i);
        srcmap = srcmap + Imap;

        if ismember(i,idx_print)
            print_count = print_count + 1;
            disp(sprintf('stack %d %% sources',print_count*5));
        end
    end
end
map = rebin_map_coarse(srcmap,10);
map = map.*100;
srcmap_out = map;

%%%%%%%%%%%%%%%%%%%%%%%%%% make ihlmap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get (stacked) psf params %%%
psfdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

bestparam = fitpsfdat(ifield).bestparam;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);
npix = 1200;
radmap = make_radius_map(zeros(2*npix+1),npix,npix).*0.7;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
psfmap = psfmap./sum(psfmap(:));
ihlamp = psfmap(npix+1,npix+1);

%%% make srcmap %%%
Nsrc = numel(subI_arr);

srcmap=zeros(7200);

disp(sprintf('making IHL map TM%d %s, mrange=(%d,%d), %d srcs',...
    inst,dt.name,m_min,m_max,Nsrc));

Nlarge = 10240+300+300;
radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;

model1 =  @(p,x) p(1).*x.^p(2);
y_r2 = model1(params(1:2),radius(2));
model2 =  @(p,x) y_r2.*(x/radius(2)).^p;

ssp1 = find(radmap<radius(2));
ssp2 = find(radmap>=radius(2) & radmap<radius(3));

Imap_large = zeros(size(radmap));
Imap_large(ssp1) = ihlamp.*model1(params(1:2),radmap(ssp1));
Imap_large(ssp2) = ihlamp.*model2(params(3),radmap(ssp2));

clear ssp1 ssp2

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
    Imap = Imap_large(1+dx:7200+dx,1+dy:7200+dy).*subI_arr(i);

    srcmap = srcmap + Imap;

    if ismember(i,idx_print)
        print_count = print_count + 1;
        disp(sprintf('make %d %% sources',print_count*5));
    end
end
map = rebin_map_coarse(srcmap,10);
ihlmap_out = map;

%%%%%%%%%%%%%%%%%%%%%%%%%% make mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rad_arr = -3.08e4.*exp(-(subm_arr-10.84).^2./(5.253).^2) + ...
    3.091e4.*exp(-(subm_arr-10.83).^2./(5.26).^2); % [arcsec]
rad_arr(find(rad_arr<7))=7;

mask = ones(720);

if numel(subm_arr)>100
    idx_print=floor(numel(subm_arr)/100) * (1:100);
else
    idx_print=[];
end
print_count=0;  

for i=1:numel(subm_arr)
    radmap=make_radius_map(mask,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7);
    mask(sp1)=0;

    if ismember(i,idx_print)
        print_count = print_count + 1;
        disp(sprintf('making mask for %d %% sources',print_count));
    end

end

end