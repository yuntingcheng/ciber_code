function [Imean_arr,Istd_arr]=...
    srcflux_2m(flight,inst,ifield,m_min,m_max,dx_arr,cbmap,mask_inst,stackband)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stack src based on PanSTARRS catalog
%
%Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - type: 1:gal, -1:star
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - dx: stamp size is 2*dx+1, dx is 10 times finer CIBER pix
% - cbmap:
% - mask_inst:
% - verbose:1 or 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PSC/');
dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.txt');

M = csvread(catfile,1);
x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

if strcmp(stackband,'y')
    m_arr=squeeze(M(:,10)');
elseif strcmp(stackband,'Ilin')
    lambdaeff=1.05;
    m_arr=squeeze(M(:,8)');
elseif strcmp(stackband,'Hlin')
     lambdaeff=1.79;
    m_arr=squeeze(M(:,9)');
end
sr = ((7./3600.0)*(pi/180.0)).^2;    
I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);


sp=find(x_arr>0.5 & x_arr<1024.5 & y_arr>0.5 & y_arr<1024.5 ...
        & m_arr<=m_max & m_arr>m_min);

x_arr = x_arr(sp);
y_arr = y_arr(sp);
m_arr = m_arr(sp);
I_arr = I_arr(sp);

%%% set up stacking %%%

idx_arr = 1:numel(m_arr);

cbmap = cbmap.*mask_inst;
dxmax = round((max(dx_arr) - 1)/2);
cbmap1 = padarray(cbmap,[dxmax dxmax],0,'both');
mask1 = padarray(mask_inst,[dxmax dxmax],0,'both');
Iint_arr = NaN([numel(m_arr), numel(dx_arr)+1]);
Iint_arr(:,end) = I_arr(:);
for i=idx_arr
    xi = round(x_arr(i)) + dxmax;
    yi = round(y_arr(i)) + dxmax;    
    for j = 1:numel(dx_arr)
        dx = round((dx_arr(j) - 1)/2);
        stamp = cbmap1(xi - dx : xi + dx, yi - dx : yi + dx);
        mstamp = mask1(xi - dx : xi + dx, yi - dx : yi + dx);
        if ~any(~mstamp(:))
            Iint_arr(i,j) = sum(stamp(:));
        end
    end
end
Imean_arr = nanmedian(Iint_arr);
Istd_arr = nanstd(Iint_arr);
disp(sprintf('count %d src between %.1f<m<%.1f '...
        ,numel(idx_arr),m_min,m_max));

return
