%dt = 6.8133e-6;
%N=6815744;
%time=(0:(N-1))*dt;
%Fs = 1./dt;
q=1;
binfac=10;
Fs = cp(inst).framerate*((1024/2)^2);
dt=1/Fs;
N = length(ts(1).t);
timefull=(0:(N-1))*dt;
lowtime = 0:dt*binfac:max(timefull(:));
sp = find(ts(q).t ~= 0);
lowy = interp1(timefull(sp),ts(q).t(sp),lowtime);
good = find((lowy == lowy));
ts(q).biny=lowy(good);
ts(q).binx= lowtime(good);

dt=ts(q).binx(2)-ts(q).binx(1);
time=ts(q).binx;
Fs=1./dt;
%% do polyfit
npoly=5;
p=polyfit(ts(q).binx,ts(q).biny,npoly);
v = polyval(p,time);
v = polyval(p,ts(q).binx);
ts(q).biny = ts(q).biny(:) - v(:);
%%
ntones = 1;

nu1 = 9.5;

nus = nu1*(1:ntones);
omegas = nus*2*pi;
phi =.5;

y=sin(time*omegas(1) + phi);
for i=1:ntones
    y=y+sin(time*omegas(i) + phi);
end

noise = awgn_pk(y*0,.04);

slope = .2;

y = y+noise +slope.*time;
%y=y+sin(time*omega1*2);

%%

dirty = y;
%dirty=ts(q).t;
dirty=ts(q).biny;

for i=1:ntones
    d = designfilt('bandstopiir','FilterOrder',2, ...
        'HalfPowerFrequency1',nus(i)*.99,'HalfPowerFrequency2',nus(i)*1.01, ...
        'DesignMethod','butter','SampleRate',Fs);
    
    clean = filtfilt(d,dirty);
    dirty=clean;
end
%%
figure
plot(time,y);
hold on
plot(time,clean);
plot(time,y);

hold off
%%
figure
[nu,bef] = philtimeps(time,y);
[nu,aft] = philtimeps(time,clean);

loglog(nu,bef);
hold on
loglog(nu,aft);
hold off
%%
lowresid = ts(q).biny-clean;
hiresid = interp1(lowtime(good),lowresid,timefull,'spline');

[nu,bef] = philtimeps(time,ts(q).t);
[nu,aft] = philtimeps(time,hiresid);
figure
loglog(nu,bef);
hold on
loglog(nu,aft);
hold off

