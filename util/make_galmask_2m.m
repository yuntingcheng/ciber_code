function mask=make_galmask_2m(flight,inst,ifield,m_min,m_max)
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

mask=zeros(1024);
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

cls_arr=squeeze(M(:,12)');


sp=find(m_arr<=m_max & m_arr>m_min & cls_arr==1);
subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);
rad_arr=1930.*subm_arr.^-0.8053-172.6; % [arcsec]

qmask=ones(512);
for i=1:numel(subm_arr)
    radmap=make_radius_map(qmask,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < rad_arr(i)./7);
    qmask(sp1)=0;
end

%%% stick to large mask %%%%
if iquad==1
    mask(1:512,1:512)=qmask;
elseif iquad==2
    mask(513:1024,1:512)=qmask;
elseif iquad==3
    mask(1:512,513:1024)=qmask;
elseif iquad==4
    mask(513:1024,513:1024)=qmask;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%%% hand mask UK no SWIRE data corners %%%
if ifield==8
    mask(1:120,800:end)=0;
    mask(980:end,800:end)=0;
end

end
