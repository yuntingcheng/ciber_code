function strmask=tm_strmask_trial(flight,inst,ifield,m_min,m_max,varargin)

  p = inputParser;
  
  p.addRequired('flight',@isnumeric);
  p.addRequired('inst',@isnumeric);
  p.addRequired('ifield',@isnumeric);
  p.addRequired('m_min',@isnumeric);
  p.addRequired('m_max',@isnumeric);
  
  p.parse(flight,inst,ifield,m_min,m_max,varargin{:});

  flight = p.Results.flight;
  inst = p.Results.inst;
  ifield = p.Results.ifield;
  m_min = p.Results.m_min;
  m_max = p.Results.m_max;
  
  clear p varargin;
  
quad_arr=['A','B','C','D'];
dt=get_dark_times(flight,inst,ifield);

catdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PSC/');

strmask=zeros(1024);

for iquad=1:4
    quad=quad_arr(iquad);
    catfile=strcat(catdir,dt.name,'_',quad,'_2m.txt');

M = csvread(catfile,1);

x_vec=squeeze(M(:,5)');
y_vec=squeeze(M(:,4)');
if inst==1
    m_vec=squeeze(M(:,6)');
else
    m_vec=squeeze(M(:,7)');
end

rad_vec=zeros(1,numel(m_vec));

for i=1:numel(rad_vec)
    if (m_vec(i)>m_min) & (m_vec(i)<m_max)
    rad_vec(i)=17;
    else
    rad_vec(i)=0;
    end
end

sp=find(m_vec<m_max);
subm_vec=m_vec(sp);
subrad_vec=rad_vec(sp);
subx_vec=x_vec(sp);
suby_vec=y_vec(sp);

mask=ones(512);

for j=1:numel(subm_vec)
    radmap=make_radius_map(mask,subx_vec(j),suby_vec(j));
    sp1 = find (radmap < subrad_vec(j));
    mask(sp1)=0;
end

if iquad==1
    strmask(1:512,1:512)=mask;
elseif iquad==2
    strmask(513:1024,1:512)=mask;
elseif iquad==3
    strmask(1:512,513:1024)=mask;
elseif iquad==4
    strmask(513:1024,513:1024)=mask;
end

end

end