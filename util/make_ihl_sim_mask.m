function [mask_all,mask_use] = make_ihl_sim_mask...
    (flight,inst,ifield,corr,m_min,m_max)
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

catdir = strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/corr_sim/TM',...
    num2str(inst),'/');

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

%%% select cat data %%%
sp0=find(m_arr<=m_max & m_arr>m_min & use_arr==0);
sp1=find(m_arr<=m_max & m_arr>m_min & use_arr==1);

subm_arr0=m_arr(sp0);
subx_arr0=x_arr(sp0);
suby_arr0=y_arr(sp0);

subm_arr1=m_arr(sp1);
subx_arr1=x_arr(sp1);
suby_arr1=y_arr(sp1);


%%% make srcmap %%%
Nsrc = numel(subm_arr0) + numel(subm_arr1);
Nsrc0 = numel(subm_arr0);
Nsrc1 = numel(subm_arr1);

rad_arr1 = -3.08e4.*exp(-(subm_arr1-10.84).^2./(5.253).^2) + ...
    3.091e4.*exp(-(subm_arr1-10.83).^2./(5.26).^2); % [arcsec]
rad_arr1(find(rad_arr1<7))=7;

rad_arr0 = -3.08e4.*exp(-(subm_arr0-10.84).^2./(5.253).^2) + ...
    3.091e4.*exp(-(subm_arr0-10.83).^2./(5.26).^2); % [arcsec]
rad_arr0(find(rad_arr0<7))=7;

mask = ones(1024);

if Nsrc>10
    idx_print=floor(Nsrc/10) * (1:10);
else
    idx_print=[];
end
print_count=0;  

for i=1:Nsrc1
    radmap=make_radius_map(mask,subx_arr1(i),suby_arr1(i));
    sp = find (radmap < rad_arr1(i)./7);
    mask(sp)=0;

    if ismember(i,idx_print)
        print_count = print_count + 1;
        disp(sprintf('making mask for %d %% sources',print_count*10));
    end
end

mask_use = mask;
for i=1:Nsrc0
    radmap=make_radius_map(mask,subx_arr0(i),suby_arr0(i));
    sp = find (radmap < rad_arr0(i)./7);
    mask(sp)=0;

    if ismember(i+Nsrc1,idx_print)
        print_count = print_count + 1;
        disp(sprintf('making mask for %d %% sources',print_count*10));
    end
end
mask_all = mask;

end