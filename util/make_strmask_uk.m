function strmask=make_strmask_uk(flight,inst,ifield,alpha,beta,m_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce a star mask from catalog
% the masking radius is defined by r= alpha *m + beta
%
% the magnitude data is in AB magnitude, and linear interpolate to 
% CIBER I and H effective wavelength. 
%
%Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - alpha: radius slope [arcsec/mag]
% - beta: radius intersection [arcsec]
% - m_max: max masking magnitude
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/UKIDSS/');
quad_arr=['A','B','C','D'];
dt=get_dark_times(flight,inst,ifield);

strmask=zeros(1024);
for iquad=1:4

%%% read cat data %%%
quad=quad_arr(iquad);
catfile=strcat(catdir,dt.name,'_',quad,'_uk.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,5)');
y_arr=squeeze(M(:,4)');
x_arr=x_arr+1;
y_arr=y_arr+1;

if inst==1
    m_arr=squeeze(M(:,6)');
else
    m_arr=squeeze(M(:,7)');
end

if ifield==8
    pn_arr=squeeze(M(:,14)');
else
    pn_arr=squeeze(M(:,15)');
end

%%% produce mask %%%%
rad_arr=alpha.*m_arr+beta; % [arcsec]
rad_arr=rad_arr./7; % [pix]

% in case r goes to negative before m_max
rad_arr(find(rad_arr<1))=1;

sp=find(m_arr<m_max & pn_arr<0.04);
subm_arr=m_arr(sp);
subrad_arr=rad_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

mask=ones(512);

for i=1:numel(subm_arr)
    radmap=make_radius_map(mask,subx_arr(i),suby_arr(i));
    sp1 = find (radmap < subrad_arr(i));
    mask(sp1)=0;
end

%%% stick to large mask %%%%
if iquad==1
    strmask(1:512,1:512)=mask;
elseif iquad==2
    strmask(513:1024,1:512)=mask;
elseif iquad==3
    strmask(1:512,513:1024)=mask;
elseif iquad==4
    strmask(513:1024,513:1024)=mask;
end

%%% hand mask UK no SWIRE data corners %%%
if ifield==8
    strmask(1:120,800:end)=0;
    strmask(980:end,800:end)=0;
end


end

end