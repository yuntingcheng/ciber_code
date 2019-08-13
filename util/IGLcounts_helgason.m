function N = IGLcounts_helgason(band,mag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the galaxy number counts from Helgason model
% data retrived from Fig.5 of
% http://adsabs.harvard.edu/abs/2012ApJ...752..113H
% 
% band: 1(2) use 1.25(1.63) data in the paper
% mag: interpolate to the given magnitude (AB)
%
% N: N/mag/deg^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

datadir = '/Users/ytcheng/ciber/doc/20170904_External/helgason/';

if band==1
    M = csvread(strcat(datadir,'Helgason125.txt'),1);
elseif band==2
    M = csvread(strcat(datadir,'Helgason163.txt'),1);
end

logN = interp1(M(:,1),log(M(:,2)),mag);
N = exp(logN);

end
