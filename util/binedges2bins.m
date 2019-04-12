function bins = binedges2bins(binedges, varargin)
  p = inputParser;
  
  p.addRequired('binedges',@isnumeric);
  p.addOptional('log',false,@islogical);
  
  p.parse(binedges,varargin{:});

  binedges = p.Results.binedges;
  log = p.Results.log;
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if log
    bins = sqrt(binedges(1:end-1) .* binedges(2:end));
else
    bins = (binedges(1:end-1) + binedges(2:end))./2;
end

return