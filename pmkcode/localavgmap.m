function [avgmap,noisemap,indexmap] = localavgmap(map,mask,nblocks)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Divide the 2D map into nblocks x nblocks chunks, and calculate the mean
%and std of non-masked region in each chunks
%
%Input:
% - map
% - mask
% - nblocks: # of blocks in both dimension
%
%Output:
% (output map has the same dimension as input, but each pixels in the same
% chuck has the same values in all three output map, so the independent
% dimension in output is only nblock x nblock)
% - avgmap: ave of each chunck
% - noisemap: std of chunk
% - indexmap: the index label of the chuncks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up resolutions
[nx,ny]= size(map);
avgmap=zeros([nx,ny]);noisemap=zeros([nx,ny]);indexmap=zeros([nx,ny]);
dx=round(nx/nblocks);dy=round(ny/nblocks);
count = 0;
%% do the loop
for i =1:(nblocks)     
        xmin=(i-1).*dx;
        if xmin == 0
            xmin =1;
        end
        xmax=(i).*dx;
         if xmax > nx
                xmax =nx;
         end
        for j=1:(nblocks)
            count = count + 1;
            ymin=(j-1).*dy;
            if ymin == 0 
                ymin =1;
            end
            ymax=(j).*dy;
            if ymax > ny
                ymax =ny;
            end
            mapchunk=map((xmin:xmax),(ymin:ymax));
            maskchunk=mask((xmin:xmax),(ymin:ymax));
            mapchunk=mapchunk(find(maskchunk));
            avgmap((xmin:xmax),(ymin:ymax))=nanmean(mapchunk);
            noisemap((xmin:xmax),(ymin:ymax))=std(mapchunk);
            indexmap((xmin:xmax),(ymin:ymax))=count;
        end
end
return
