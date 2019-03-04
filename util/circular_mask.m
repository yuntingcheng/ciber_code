function mask=circular_mask(x,y,r,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%
%generate a circle mask centered at [x,y] w/ radius r.
%Input:
%x,y: circle center coordinate
%r: circle radius
%maskin(Optional): input mask, defulat: ones(1024)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse data
  p = inputParser;
  
  p.addRequired('x',@isnumeric);
  p.addRequired('y',@isnumeric);
  p.addRequired('r',@isnumeric);
  p.addOptional('maskin',ones(1024),@isnumeric);
  
  p.parse(x,y,r,varargin{:});

  x=p.Results.x;
  y=p.Results.y;
  r=p.Results.r;
  maskin=p.Results.maskin;
  clear p varargin;
%%
[Nx,Ny]=size(maskin);
[xx,yy]=meshgrid(1:Nx,1:Ny);
rmap=sqrt((xx-x).^2+(yy-y).^2);
mask=maskin;mask(find(rmap<r))=0;
return
