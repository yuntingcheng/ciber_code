function [filtts,stack,stackfit,sinparams,tchu] = pulsar_filtsin(time,meas,freq,varargin)
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('time');
p.addRequired('meas');
p.addRequired('freq');
p.addOptional('makeplot',0,@isnumeric); 

p.parse(time,meas,freq,varargin{:});

time = p.Results.time;
meas = p.Results.meas;
freq = p.Results.freq;
makeplot = p.Results.makeplot;


time = time(:);
meas = meas(:);
mask = meas*0 +1;
sp = find(meas == 0);
mask(sp) = 0 ;


dt = time(2) -time(1);
period = 2./freq;
nsamp = round(period/dt);
r = nsamp - (period/dt);
stack = (1:nsamp)*0;
nchunk = round(numel(time)/nsamp) - 1;
stack = stack(:);
maskstack = stack(:);
for i=1:nchunk
    rs =i*r;
    c=round(rs);
    %stack = stack + meas(((i-1)*nsamp+1 -c): ((i)*nsamp-c));
    zz = meas(((i-1)*nsamp+1 -c): ((i)*nsamp-c));
    stack = stack+zz;
    zzm = mask(((i-1)*nsamp+1 -c): ((i)*nsamp-c));
    maskstack = maskstack+zzm;
    
end
stack = stack./maskstack;

%Starting = [(max(stack(:)) -min(stack(:))),freq,dt,.5];
%options=optimset('Display','iter');
%Estimates=fminsearch(@mysinfit,Starting,options,stack,stack);
tchu =time(1:nsamp);

phases = 0:.001:2*pi;
chi2 = phases*0;
for p=1:numel(phases)
    mo=sin(2*pi*freq.*tchu +phases(p));
    diff = (stack(:) - mo(:)).^2;
    chi2(p) = sum(diff(:));
end
[mm,sp]=min(chi2(:));
bestphase = phases(sp);
mo=sin(2*pi*freq.*tchu +bestphase);
line = linfit(mo,stack);
bestamp = line(1);
bestoff = line(2);

sinparams = struct('amp',bestamp,'off',bestoff,'phase',bestphase,'freq',freq);

stackfit = bestamp*sin(2*pi*freq.*tchu +bestphase)+bestoff;


if makeplot
    plot(tchu,stack,'color','black')
    hold on
    plot(tchu,stackfit,'linewidth',2,'color','red')
    hold off
end

bigmo = bestamp*sin(2*pi*freq.*time +bestphase)+bestoff;
filtts = meas(:)-bigmo(:);

return