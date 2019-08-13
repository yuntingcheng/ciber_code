function [mask,num]=make_mask_sim(flight,inst,ifield,m_min,m_max)
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
dt=get_dark_times(flight,inst,ifield);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');
catfile=strcat(catdir,'sides.txt');
M = csvread(catfile,1);
x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
xg_arr=x_arr+1;
yg_arr=y_arr+1;

if inst==1
    mg_arr = squeeze(M(:,6)');
elseif inst==2
    mg_arr = squeeze(M(:,7)');
end

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/trilegal/');
catfile=strcat(catdir,dt.name, '.csv');
M = csvread(catfile,1);
x_arr=squeeze(M(:,3)');
y_arr=squeeze(M(:,4)');
xs_arr=x_arr+1;
ys_arr=y_arr+1;

if inst==1
    ms_arr = squeeze(M(:,1)');
elseif inst==2
    ms_arr = squeeze(M(:,2)');
end

x_arr = [xs_arr, xg_arr];
y_arr = [ys_arr, yg_arr];
m_arr = [ms_arr, mg_arr];


sp=find(m_arr<=m_max & m_arr>m_min);
subm_arr=m_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

rad_arr = get_mask_radius(inst,ifield,subm_arr);

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
        fprintf('making mask for %d %% sources\n',print_count);
    end

end

return