function mask = stamp_clip2(map,max_arr,min_arr,xcent,ycent,nbins,rmin)
mask = zeros(size(map));
mask(find(map)) = 1;
rad = make_radius_map(map,xcent,ycent);
binedge=logspace(log10(rmin),log10(max(rad(:))),nbins+1);
if rmin==0 ; binedge(1)=0 ; end

for ibin=1:nbins
    clipmax = max_arr(ibin);
    clipmin = min_arr(ibin);

    spclip = find(rad >= binedge(ibin) ...
        & rad < binedge(ibin+1) & mask ...
        & (map >= clipmax | map <= clipmin));
    
    clip_val = unique(map(spclip))';
    
    mask(ismember(map,clip_val)) = 0;
    
end

return