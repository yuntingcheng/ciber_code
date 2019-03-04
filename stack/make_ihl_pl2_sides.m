function [map] = make_ihl_pl2_sides(flight,inst,ifield,interp,...
    m_min,m_max,params,radius)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the sim IHL power-law profile
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - rcut: cutoff radius if slope > -2
% - m_min: min masking magnitude (PS y band)
% - m_max: max masking magnitude (PS y band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');
psfdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');

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
sp=find(m_arr<=m_max & m_arr>m_min);
subI_arr=I_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

%%% get x,y coord in small grid %%%
xsmall_arr = subx_arr.*10 - 4.5;
ysmall_arr = suby_arr.*10 - 4.5;

%%% get (stacked) psf params %%%
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

end