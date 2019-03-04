function [fitstart,fitend]=reset_finder(ts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function find the reset points
%for a given stream of data.
%The alogrithm is that if it is a reset,
%the different with the previous one sould be
%very negative, while the different between the
%next one should be very positive.
%Input:
%   -ts: 1D arr of data 
%Output:
%   -rs_arr:reset point arr
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the reset points

nfr=numel(ts);
%res_arr record the reset positions
res_arr=[];
for i=1:nfr-1
    diff(i)=ts(i+1)-ts(i);
end
med_diff=median(abs(diff));

for i=1:nfr-2
    if abs(diff(i))/med_diff>100 && ...
       abs(diff(i+1))/med_diff>100 && ...
       diff(i)<0 && diff(i+1)>0
   
            res_arr=[res_arr,i+1];
    end

end
nres=numel(res_arr);

if nres~=0
%% find the duration between the reset points

%dur_arr is the duration(time-step) between 
%each resets. One more dimension than res_arr
%because the 1st and end element is the duration 
%to the begining and the end of the data
dur_arr=zeros(1,nres+1);

for i=1:nres
    if i==1
        dur_arr(i)=res_arr(i)-1;
    else
        dur_arr(i)=res_arr(i)-res_arr(i-1);
    end
end
dur_arr(nres+1)=nfr-res_arr(nres);

%% find the fitstart and fitend
%which is the boundary of the longest duration

maxdur=find(dur_arr==max(dur_arr));
%if more than one longest duration
%choose the first one
maxdur=maxdur(1);
ndur=numel(dur_arr);
if maxdur==1
    fitstart=3;%incase the 1st reset
    fitend=res_arr(1)-1;
elseif maxdur==ndur
    fitstart=res_arr(nres)+1;
    fitend=nfr-1;%incase the last reset
else
    fitstart=res_arr(maxdur-1)+1;
    fitend=res_arr(maxdur)-1;
end


else
    fitstart=3;
    fitend=nfr-1;
end
end


