function Q = quartile_from_hist(bins, counts, quartile)
  %%%%%%%%%%%%%%%%%%%%%%
  % Calculate the quartile value from histogram
  % bins: bin value in hist
  % counts: bin counts in hist
  % quartile: [0,1] quartile value, 0.5 - median, 0.25 - 1st quartile, 0.75
  % - 3rd quartile
  %%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('bins',@isnumeric);
  p.addRequired('counts',@isnumeric);
  p.addRequired('quartile',@isnumeric);
  p.parse(bins,counts,quartile);

  bins = p.Results.bins;
  counts = p.Results.counts;
  quartile = p.Results.quartile;
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totc = sum(counts);

c = 0;
i = 0;
while c <= totc * quartile
    i = i + 1;
    c = c + counts(i);
end

if c==totc*quartile && i < numel(bins)
    Q = (bins(i) + bins(i+1))/2;
else
    Q = bins(i);
end


return