function r_arr = get_mask_radius(inst,ifield,m_arr,varargin)
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
  p.addOptional('rmin',0.8,@isnumeric);
  
  p.parse(inst,ifield,m_arr,varargin{:});

  inst     = p.Results.inst;
  ifield   = p.Results.ifield;
  m_arr    = p.Results.m_arr;
  rmin     = p.Results.rmin;
  
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (inst ==1) && (ifield == 4)
%     p = [-1.11633722e-04 2.16652113e-03 1.15547601e-01 -4.03349233e+00 3.62368050e+01];
% elseif (inst ==1) && (ifield == 5)
%     p = [2.12352147e-03 -1.31503533e-01 3.08694967e+00 -3.36491673e+01 1.51490226e+02];
% elseif (inst ==1) && (ifield == 6)
%     p = [-3.64670160e-04 9.78870631e-03 1.97962759e-01 -8.62383830e+00 7.46680846e+01];
% elseif (inst ==1) && (ifield == 7)
%     p = [5.99100976e-04 -4.61153118e-02 1.37516259e+00 -1.92516245e+01 1.09116967e+02];
% elseif (inst ==1) && (ifield == 8)
%     p = [1.29450013e-04 -1.13807961e-02 3.75402639e-01 -5.92575257e+00 3.96534566e+01];

if (inst ==1) && (ifield == 4)
    p = [-1.37117784e-03 4.44237745e-02 -3.59547593e-01 -2.00388819e+00 3.35764491e+01];
elseif (inst ==1) && (ifield == 5)
    p = [2.54298277e-03 -1.51824969e-01 3.44127294e+00 -3.63803776e+01 1.59400972e+02];
elseif (inst ==1) && (ifield == 6)
    p = [-6.99300699e-03 2.72286972e-01 -3.54707848e+00 1.38693339e+01 2.82385436e+01];
elseif (inst ==1) && (ifield == 7)
    p = [-6.99300699e-03 2.72286972e-01 -3.54707848e+00 1.38693339e+01 2.82385436e+01];
elseif (inst ==1) && (ifield == 8)
    p = [7.67575862e-16 3.88500388e-04 5.22144522e-02 -2.66817405e+00 2.96807273e+01];
    
elseif (inst ==2) && (ifield == 4)
    p = [-3.85643768e-04 7.82904459e-03 3.26402601e-01 -1.08280100e+01 8.64444140e+01];
elseif (inst ==2) && (ifield == 5)
    p = [7.12118588e-04 -5.02993629e-02 1.39138714e+00 -1.84178273e+01 1.01306248e+02];
elseif (inst ==2) && (ifield == 6)
    p = [-2.65665707e-04 3.13656931e-03 3.48836213e-01 -1.00372716e+01 7.96480138e+01];
elseif (inst ==2) && (ifield == 7)
    p = [9.14156148e-04 -6.40994218e-02 1.75385325e+00 -2.28062250e+01 1.22170242e+02];
elseif (inst ==2) && (ifield == 8)
    p = [-1.17271789e-04 1.95533168e-03 1.31302787e-01 -4.27016665e+00 3.72392345e+01];
end

if numel(m_arr)==0
    r_arr=[];
    return
end

% set the radius to rmin if the function goes below rmin
mbins = min(m_arr):0.1:max(m_arr);
r_arr = polyval(p,mbins);
mcut = find(r_arr < rmin);
if numel(mcut) > 0
    mcut = mbins(mcut(1));
else
    mcut = max(mbins) + 1;
end

% set the radius to min point if the function starts increase with mag
r_arr = polyval(p,m_arr);
r_arr(m_arr > mcut) = rmin;

sp = find(r_arr == min(r_arr));
if numel(sp)>=2
    sp = sp(1);
end
r_arr(m_arr > m_arr(sp)) = r_arr(sp);

r_arr = r_arr.*7;
return