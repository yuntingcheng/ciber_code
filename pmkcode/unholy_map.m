function fillmap=unholy_map(map,mask,initialres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% phil korngut 9/7/2011
% This program fill in the blanks of a masked image with a local mean.
% First filled the masked region with local mean of initialres, and if
% there is still zeros in the filled map, make the block size in
% localavgmap.m twice larger and fill again, and do it iteratively unitll
% all the pixels in non-zero.

%Input:
% - map
% - mask
% - initialres: # of blocks in localavgmap in the first iteration
%
%Output:
% - fillmap: the filled map with the holes filled with local mean.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make the first high resolution avg map
avgmap = localavgmap(map,mask,initialres);
%% fill in the blanks
fillmap=map;
death = find(mask == 0);
fillmap(death)=avgmap(death);
%% see how we did with the first pass
zerospot=find(fillmap ==0 | fillmap ~= fillmap);
nzeros=length(zerospot);
%% keep going till you got all the blanks
count = 0;
while nzeros > 0
    count = count +1;
    newres=round(initialres./2);%twice the block size
    death=find(fillmap ==0 | fillmap ~= fillmap);
    avgmap = localavgmap(map,mask,newres);
    fillmap(death)=avgmap(death);
    death2=find(fillmap == 0 | fillmap ~= fillmap);
    nzeros=length(death2);
    %display(strcat('pass through=',num2str(count)))
    %display(strcat('holy pixels left=',num2str(nzeros)))
    initialres=newres;
    % make an out at a cap to avoid infinite loops
    if count > 100
        nzeros = 0;
    end
end
return