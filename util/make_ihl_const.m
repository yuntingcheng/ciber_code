function ihlmap=make_ihl_const(flight,inst,ifield,m_min,m_max,amp,rad)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce a constant ihl around sources between m_min,m_max,
%out to radius r arcsec
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - amp: ihl level compared to center pixel (0.7 arcsec) of PSF
% - m_min: min masking magnitude
% - m_max: max masking magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PSC/');

quad_arr=['A','B','C','D'];
dt=get_dark_times(flight,inst,ifield);

ihlmap=zeros(1024);
for iquad=1:4

%%% read cat data %%%
quad=quad_arr(iquad);
catfile=strcat(catdir,dt.name,'_',quad,'_2m_matched.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');
x_arr=x_arr+1;
y_arr=y_arr+1;


if inst==1
    m_arr=squeeze(M(:,6)');
    lambdaeff=1.05;
else
    m_arr=squeeze(M(:,7)');
    lambdaeff=1.79;
end
sr = ((7/3600.0)*(pi/180.0))^2;
I_arr=3631.*10.^(-m_arr./2.5).*(3./lambdaeff).*1e6./(sr.*1e9);

cls_arr=squeeze(M(:,12)');

loaddir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

bestparam=fitpsfdat(ifield).bestparam_norm;
A=bestparam(1);
B=bestparam(2);
sig=bestparam(3);
r0=bestparam(4);
alpha=bestparam(5);

radmap = make_radius_map(zeros(2*800+1),800,800).*0.7;
psfmap = A*exp(-radmap.^2./2./sig^2)+B./(1+(radmap./r0).^alpha);
ihlamp = psfmap(801,801)*amp*100;


sp=find(m_arr<=m_max & m_arr>m_min & cls_arr==1);
subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);
subI_arr=I_arr(sp);
map=zeros(512);
for i=1:numel(subm_arr)
    radmap=make_radius_map(map,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad./7);
    map(sp1)=ihlamp*subI_arr(i);
end

%%% stick to large mask %%%%
if iquad==1
    ihlmap(1:512,1:512)=map;
elseif iquad==2
    ihlmap(513:1024,1:512)=map;
elseif iquad==3
    ihlmap(1:512,513:1024)=map;
elseif iquad==4
    ihlmap(513:1024,513:1024)=map;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

end
