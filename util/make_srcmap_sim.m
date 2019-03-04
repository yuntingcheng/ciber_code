function [map_all,map_use] = make_srcmap_sim(flight,inst,ifield,corr,...
    m_min,m_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the sim IHL power-law profile
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - type: 1:gal, -1:star, 0:all, 2:undefined
% - rcut: cutoff radius if slope > -2
% - m_min: min masking magnitude (PS y band)
% - m_max: max masking magnitude (PS y band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catdir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');
psfdir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/forward_model/'));

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'_simcoord.csv');

M = csvread(catfile,1);

if corr==1
    x_arr=squeeze(M(:,1)');
    y_arr=squeeze(M(:,2)');
else
    x_arr=squeeze(M(:,3)');
    y_arr=squeeze(M(:,4)');    
end

m_arr=squeeze(M(:,5)');
use_arr=squeeze(M(:,6)');

lambdaeff=0.9633;
sr = ((7./3600.0)*(pi/180.0)).^2;
I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);


%%% select cat data %%%
sp0=find(m_arr<=m_max & m_arr>m_min & use_arr==0);
sp1=find(m_arr<=m_max & m_arr>m_min & use_arr==1);

subI_arr0=I_arr(sp0);
subx_arr0=x_arr(sp0);
suby_arr0=y_arr(sp0);

subI_arr1=I_arr(sp1);
subx_arr1=x_arr(sp1);
suby_arr1=y_arr(sp1);

%%% get x,y coord in small grid %%%
xsmall_arr0 = subx_arr0.*10 - 4.5;
ysmall_arr0 = suby_arr0.*10 - 4.5;
xsmall_arr1 = subx_arr1.*10 - 4.5;
ysmall_arr1 = suby_arr1.*10 - 4.5;

%%% get psf params %%%
psfparfile = strcat(psfdir,'TM',num2str(inst),'_',dt.name,'_bestparam.txt');
params = csvread(psfparfile);
A=params(1);
B=params(2);
sig=params(3);
r0=params(4);
alpha=params(5);

%%% make srcmap %%%
Nsrc = numel(subI_arr0) + numel(subI_arr1);
Nsrc0 = numel(subI_arr0);
Nsrc1 = numel(subI_arr1);

srcmap=zeros(10240);

disp(sprintf('making sim src map TM%d %s, mrange=(%d,%d), %d srcs',...
    inst,dt.name,m_min,m_max,Nsrc));

Nlarge = 10240+300+300;
radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
Imap_large = A*exp(-radmap.^2./2./(sig)^2)+B./(1+(radmap./(r0)).^alpha);

if Nsrc>20
    idx_print=floor(Nsrc/20) * (1:20);
else
    idx_print=[];
end

print_count=0;

for i=1:Nsrc1
    xi = round(xsmall_arr1(i));
    yi = round(ysmall_arr1(i));
    dx = Nlarge + 1 - xi;
    dy = Nlarge + 1 - yi;
    Imap = Imap_large(1+dx:10240+dx,1+dy:10240+dy).*subI_arr1(i);

    srcmap = srcmap + Imap;

    if ismember(i,idx_print)
        print_count = print_count + 1;
        disp(sprintf('make %d %% sources',print_count*5));
    end
end
map_use = rebin_map_coarse(srcmap,10).*100;

for i=1:Nsrc0
    xi = round(xsmall_arr0(i));
    yi = round(ysmall_arr0(i));
    dx = Nlarge + 1 - xi;
    dy = Nlarge + 1 - yi;
    Imap = Imap_large(1+dx:10240+dx,1+dy:10240+dy).*subI_arr0(i);

    srcmap = srcmap + Imap;

    if ismember(i+Nsrc1,idx_print)
        print_count = print_count + 1;
        disp(sprintf('make %d %% sources',print_count*5));
    end
end

map_all = rebin_map_coarse(srcmap,10).*100;
end