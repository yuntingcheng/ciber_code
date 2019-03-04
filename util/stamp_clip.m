function mask = stamp_clip(map,xcent,ycent,nbins,sig,iter_clip,rmin)
mask = zeros(size(map));
mask(find(map)) = 1;
rad = make_radius_map(map,xcent,ycent);
binedge=logspace(log10(rmin),log10(max(rad(:))),nbins+1);
if rmin==0 ; binedge(1)=0 ; end

for ibin=1:nbins
    for j=1:iter_clip
        sp=find(rad >= binedge(ibin) ...
            & rad < binedge(ibin+1) & mask);
        b = map(sp);
        
        clipmax = nanmedian(b(:))+sig*std(b(:));
        clipmin = nanmedian(b(:))-sig*std(b(:));
        spclip = find(rad >= binedge(ibin) ...
            & rad < binedge(ibin+1) & mask ...
            & (map >= clipmax | map <= clipmin));
        
        clip_val = unique(map(spclip))';
        
        mask(ismember(map,clip_val)) = 0;
        
    end
end

return