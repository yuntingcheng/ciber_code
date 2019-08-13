function [mask,num]=make_mask_hsc(flight,inst,ifield,field,m_min,m_max,Ith)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the mask from sim
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
catdir=mypaths.hsccatdir;
catfile=strcat(catdir,field,'.csv');
M = csvread(catfile,1);
x_arr=squeeze(M(:,3)');
y_arr=squeeze(M(:,4)');
x_arr=x_arr-1.5;
y_arr=y_arr-1.5;
Npix_x = 1024;
Npix_y = 1024;

if inst==1
    m_arr = squeeze(M(:,10)');
elseif inst==2
    m_arr = squeeze(M(:,11)');
end


sp=find(m_arr<=m_max & m_arr>m_min);
subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

rad_arr = get_mask_radius_th(inst,ifield,subm_arr,Ith);

mask = ones(Npix_x);
num = zeros(Npix_y);

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
        fprintf('making mask for %d %% sources\n',print_count);
    end

end

return