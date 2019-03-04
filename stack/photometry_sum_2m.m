function [I1, I3, I5, IT]=photometry_sum_2m...
    (flight,inst,ifield,m_min,m_max,cbmap,mask_inst,strmask,strnum)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stack src based on 2MASS catalog
% 
% Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - cbmap:
% - mask_inst:
% - strmask:
% - strnum:
% 
% Output:
% I1, I3, I5: sum of signal in 1x1, 3x3, 5x5 pixel around the source [map unit]
% IT: true brightness from catalog mag [nW/m2/sr]
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

rad_arr = get_mask_radius(inst,ifield,M(:,10)');

if inst == 1
    m_arr = squeeze(M(:,8)');
    lambdaeff=1.05;
else
    m_arr = squeeze(M(:,9)');
    lambdaeff=1.79;
end

sp=find(x_arr>3 & x_arr<1022 & y_arr>3 & y_arr<1022);

x_arr = x_arr(sp);
y_arr = y_arr(sp);
m_arr = m_arr(sp);
rad_arr = rad_arr(sp);

%%% count the center pix map
xround_arr=round(x_arr);
yround_arr=round(y_arr);

centnum_map = zeros(1024);
for i=1:numel(xround_arr)
    centnum_map(xround_arr(i),yround_arr(i))=...
        centnum_map(xround_arr(i),yround_arr(i))+1;
end

%%% select sources %%%
sp=find(m_arr<=m_max & m_arr>m_min);

m_arr=m_arr(sp);
x_arr=x_arr(sp);
y_arr=y_arr(sp);
rad_arr = rad_arr(sp);

%%% select sources not coexist with others in the same pixel %%%
sp_use = [];
for i=1:numel(sp)
    if centnum_map(round(x_arr(i)),round(y_arr(i)))==1 ...        
            & mask_inst(round(x_arr(i)),round(y_arr(i)))==1
        sp_use = [sp_use, i];
    end
end

m_arr = m_arr(sp_use);
x_arr = x_arr(sp_use);
y_arr = y_arr(sp_use);
rad_arr = rad_arr(sp_use);

sr = ((7./3600.0)*(pi/180.0)).^2;
I_arr = 3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);

%%% set up stacking %%%
I1 = [];
I3 = [];
I5 = [];
IT = [];

idx_arr = 1:numel(m_arr);


for i=idx_arr
    x = round(x_arr(i));
    y = round(y_arr(i));
    cbmapi = cbmap.*strmask.*mask_inst;
    
    %%% unmask the target
    radmap = make_radius_map(cbmapi,x_arr(i),y_arr(i));
    %sp1 = find(radmap < rad_arr(i)./7 & strnum==1 & mask_inst==1);
    sp1 = find(radmap < rad_arr(i)./7 & mask_inst==1);
    cbmapi(sp1) = cbmap(sp1);
    
    bigstamp = cbmapi(x-2: x+2, y-2: y+2);
    if ~all(bigstamp(:))
        continue
    end
    
    I1 = [I1, cbmapi(x,y)];
    I3 = [I3, sum(sum(cbmapi(x-1: x+1, y-1: y+1)))];
    I5 = [I5, sum(sum(cbmapi(x-2: x+2, y-2: y+2)))];
    IT = [IT, I_arr(i)];

end
%%

%disp(sprintf('avg %d src between %.1f<m<%.1f ',numel(IT),m_min,m_max));

return
