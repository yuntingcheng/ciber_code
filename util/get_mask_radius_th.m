function r_arr = get_mask_radius_th(inst,ifield,m_arr,Ith,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masking function
% Input:
% m_arr - I/H band magnitude from catalog.
% r_arr - masking radius in arcsec
%
% r_arr is set to 3.5 if it is <3.5 from the fitting funciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  p = inputParser;
  
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('m_arr',@isnumeric);
  p.addRequired('Ith',@isnumeric);
  p.addOptional('rmin',0.8,@isnumeric);
  
  p.parse(inst,ifield,m_arr,Ith,varargin{:});

  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  m_arr    = p.Results.m_arr;
  Ith      = p.Results.Ith;
  rmin     = p.Results.rmin;
  
  clear p varargin;
  
flight=40030;
mypaths=get_paths(flight);

loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

beta=fitpsfdat(ifield).psfmodel.beta_best;
rc=fitpsfdat(ifield).psfmodel.rc_best;
norm=fitpsfdat(ifield).psfmodel.norm;

Nlarge = 100;
radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
Imap_large = norm .* (1 + (radmap/rc).^2).^(-3.*beta./2);

if inst==1
    lambdaeff=1.05;
else
    lambdaeff=1.79;
end
sr = ((7./3600.0)*(pi/180.0)).^2;
I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
r_arr = zeros(size(m_arr));
for im=1:numel(m_arr)
    I = I_arr(im);
    sp = find(Imap_large.*I > Ith./100);
    if numel(sp)>0
        r_arr(im) = max(radmap(sp));
    end
end


if numel(m_arr)==0
    r_arr=[];
    return
end

if ~isnan(rmin)
    r_arr(r_arr<rmin)=rmin.*7;
end

return