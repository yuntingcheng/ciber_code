function [res_arr]=get_resets(ts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function find the reset points for a given stream of data.
%The alogrithm is that if it is a reset, the different with the 
%previous one sould be very negative, while the different between
%the next one should be very positive.
%Input:
%   -ts: 1D arr of data 
%Output:
%   -res_arr:reset point arr 1->reset, 0->not reset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the reset points

nfr=numel(ts);
%res_arr record the reset positions
res_arr=zeros(1,nfr);
for i=1:nfr-1
    diff(i)=ts(i+1)-ts(i);
end
med_diff=median(abs(diff));

for i=1:nfr-2
    if abs(diff(i))/med_diff>100 && ...
       abs(diff(i+1))/med_diff>100 && ...
       diff(i)<0 && diff(i+1)>0   
            
       res_arr(i+1)=1;
    end

end
end


