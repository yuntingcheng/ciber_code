function [fit_map,mx,my]=plane_fit(map,mask,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Given a map and mask, fit a plane with least square. 
%Input:
%   - map
%   - mask
%   - (optional)iter: # of iteration, default:10
%Output:
%   - fit_map: the plane fitted
%   -mx:grad in x direction across the map
%   -my:grad in y direction across the map
%
%Algorithm:
%here we fit a plane with #iter of iterrations. 
%We first take the 2*(max-min) of the masked map be the limit of slope, 
%and divide into 11 numbers in x and y respectively, then find out the min
%chi2.  Then take the region between best fit slope +- 1 grid in first 
%iter as the limit of 2nd iter, and again divide into 11 grids and find
%chi2. Then do the 3rd iteration. 
%
%NOTE: MUST MASK OUT THE CRAZY PIXEL! Since we take the min and max of 
%masked map as the 1st iter limit, this breaks if there is crazy pixels 
%not masked. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%parse data
p=inputParser;
p.addRequired('map',@isnumeric);
p.addRequired('mask',@isnumeric);
p.addOptional('iter',10,@isscalar);
p.parse(map,mask,varargin{:});
map=p.Results.map;
mask=p.Results.mask;
iter=p.Results.iter;
clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%

[nx ny]=size(map);
[xx,yy]=meshgrid(1:nx,1:ny);
xx=xx./nx;xx=xx-mean(xx(:));
yy=yy./ny;yy=yy-mean(yy(:));

map(find(map~=map))=0;
map(find(map==inf))=0;
map(find(map==-inf))=0;
mm=map;mm(find(mask==0))=0;
mm=mm-mean(mm(find(mm)));mm=mm.*mask;

% start 1st iteration
%pr=sprintf('plane_fit iter #%d',1);disp(pr);
max1=max(map(find(map)));min1=min(map(find(map)));
m1_arr=linspace((min1-max1)*2,(max1-min1)*2,11);
mx_b=nan;my_b=nan;chi2_b=inf;plane_b=zeros(size(map));i_b=0;j_b=0;
for i=1:11
    for j=1:11
        mx=m1_arr(i);my=m1_arr(j);
        plane=mx.*xx+my.*yy;mp=plane.*mask;
        chi2_map=(mp-mm).^2;chi2=sum(chi2_map(:));
        if chi2<chi2_b
            chi2_b=chi2;mx_b=mx;my_b=my;plane_b=plane;i_b=i;j_b=j;
        end
    end
end

if i_b==1|i_b==11|j_b==1|j_b==11
    warning('plane_fit iter 1 out of range')
end

% start 2nd iteration
%pr=sprintf('plane_fit iter #%d',2);disp(pr);
max2x=m1_arr(i_b+1);min2x=m1_arr(i_b-1);
max2y=m1_arr(j_b+1);min2y=m1_arr(j_b-1);
mx2_arr=linspace(min2x,max2x,11);my2_arr=linspace(min2y,max2y,11);
mx_b=nan;my_b=nan;chi2_b=inf;plane_b=zeros(size(map));i_b=0;j_b=0;
for i=1:11
    for j=1:11
        mx=mx2_arr(i);my=my2_arr(j);
        plane=mx.*xx+my.*yy;mp=plane.*mask;
        chi2_map=(mp-mm).^2;chi2=sum(chi2_map(:));
        if chi2<chi2_b
            chi2_b=chi2;mx_b=mx;my_b=my;plane_b=plane;i_b=i;j_b=j;
        end
    end
end

if i_b==1|i_b==11|j_b==1|j_b==11
     warning('plane_fit iter 2 out of range')
end

% start >3rd iterations
for iteration=3:iter
%pr=sprintf('plane_fit iter #%d',iteration);disp(pr);

max2x=mx2_arr(i_b+1);min2x=mx2_arr(i_b-1);
max2y=my2_arr(j_b+1);min2y=my2_arr(j_b-1);
mx2_arr=linspace(min2x,max2x,11);my2_arr=linspace(min2y,max2y,11);
mx_b=nan;my_b=nan;chi2_b=inf;plane_b=zeros(size(map));i_b=0;j_b=0;
for i=1:11
    for j=1:11
        mx=mx2_arr(i);my=my2_arr(j);
        plane=mx.*xx+my.*yy;mp=plane.*mask;
        chi2_map=(mp-mm).^2;chi2=sum(chi2_map(:));
        if chi2<chi2_b
            chi2_b=chi2;mx_b=mx;my_b=my;plane_b=plane;i_b=i;j_b=j;
        end
    end
end

if i_b==1|i_b==11|j_b==1|j_b==11
    warning('plane_fit iter >3 out of range')
end
end

fit_map=plane_b;mx=mx_b;my=my_b;

return