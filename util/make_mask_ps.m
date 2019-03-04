function [mask,num]=make_mask_ps(flight,inst,band,ifield,type,m_min,m_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the mask from PanSTARR
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

catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PanSTARRS/');

dt=get_dark_times(flight,inst,ifield);
%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

cls_arr=squeeze(M(:,10)');

if band == 'y'
    m_arr=squeeze(M(:,9)');
elseif band == 'Ilin'
    m_arr=squeeze(M(:,11)');
elseif band == 'Hlin'
    m_arr=squeeze(M(:,12)');
elseif band == 'I'
    my_arr=squeeze(M(:,9)');
    mlin_arr=squeeze(M(:,11)');
    [m_arr, ~] = get_corrected_mag(inst, mlin_arr, my_arr, cls_arr);
elseif band == 'H'
    my_arr=squeeze(M(:,9)');
    mlin_arr=squeeze(M(:,12)');
    [m_arr, ~] = get_corrected_mag(inst, mlin_arr, my_arr, cls_arr);
end

sp=find(m_arr<=m_max & m_arr>m_min & cls_arr==type);

if type == 0
    sp=find(m_arr<=m_max & m_arr>m_min);
end

if type == 2
    sp=find(m_arr<=m_max & m_arr>m_min & cls_arr~=1 & cls_arr~=-1);
end

subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

rad_arr = get_mask_radius(inst,ifield,subm_arr);

mask = ones(1024);
num = zeros(1024);

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
    num(sp1) = num(sp1) +1;

    if ismember(i,idx_print)
        print_count = print_count + 1;
        disp(sprintf('making mask for %d %% sources',print_count));
    end

end

end
