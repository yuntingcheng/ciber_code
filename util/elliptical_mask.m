function mask=elliptical_mask(x,y,a,b,theta,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%
%generate a elliptical mask centered at [x,y],
%semi-major/minor axis [a,b],
%inclination angle (a wrt x axis) theta degree (theta=[0,180])
%
%Input:
%x,y: circle center coordinate
%a,b: semi-major/minor axis
%theta: inclination angle in degree
%maskin(Optional): input mask, defulat: ones(1024)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parse data
  p = inputParser;
  
  p.addRequired('x',@isnumeric);
  p.addRequired('y',@isnumeric);
  p.addRequired('a',@isnumeric);
  p.addRequired('b',@isnumeric);
  p.addRequired('theta',@isnumeric);
  p.addOptional('maskin',ones(1024),@isnumeric);
  
  p.parse(x,y,a,b,theta,varargin{:});

  x=p.Results.x;
  y=p.Results.y;
  a=p.Results.a;
  b=p.Results.b;
  theta=p.Results.theta;
  maskin=p.Results.maskin;
  clear p varargin;
%%
[Nx,Ny]=size(maskin);
[xx,yy]=meshgrid(1:Nx,1:Ny);
xxnew=xx.*cos(theta*pi/180)+yy.*sin(theta*pi/180);
yynew=-xx.*sin(theta*pi/180)+yy.*cos(theta*pi/180);

xnew=x*cos(theta*pi/180)+y*sin(theta*pi/180);
ynew=-x*sin(theta*pi/180)+y*cos(theta*pi/180);

rmap=((xxnew-xnew)/a).^2+((yynew-ynew)/b).^2;
mask=maskin;mask(find(rmap<1))=0;
return
