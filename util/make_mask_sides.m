function [mask,num]=make_mask_sides(flight,inst,m_min,m_max)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the mask from SIDES catalog
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - m_min: min masking magnitude (PS y band)
% - m_max: max masking magnitude (PS y band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mypaths=get_paths(flight);
catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');

%%% read cat data %%%
catfile=strcat(catdir,'sides.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

m_arr=squeeze(M(:,4)');

%%% select cat data %%%
sp=find(m_arr<=m_max & m_arr>m_min);
subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

rad_arr = -3.08e4.*exp(-(subm_arr-10.84).^2./(5.253).^2) + ...
    3.091e4.*exp(-(subm_arr-10.83).^2./(5.26).^2); % [arcsec]
rad_arr(find(rad_arr<7))=7;

mask = ones(720);
num = zeros(720);

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
