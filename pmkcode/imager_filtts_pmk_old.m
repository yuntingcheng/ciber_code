function [filtmap,dumbmap,sigmask] = imager_filtts(fr1,varargin)
 
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('fr1');
p.addOptional('domapin',1,@isnumeric); 
p.addOptional('domaskin',1,@isnumeric); 
p.addOptional('sigma',3,@isnumeric); 
p.addOptional('npoly',12,@isnumeric); 
p.addOptional('binfac',100,@isnumeric); 
p.addOptional('dopoly',1,@isnumeric); 
p.addOptional('donotch',0,@isnumeric); 
p.addOptional('filtwidth',.01,@isnumeric); 
p.addOptional('mapin',ones(1024),@isnumeric); 
p.addOptional('maskin',ones(1024),@isnumeric);
p.addOptional('offin',ones(1024),@isnumeric); 

p.parse(fr1,varargin{:});

domaskin = p.Results.domaskin;
domapin = p.Results.domapin;
maskin = p.Results.maskin;
mapin = p.Results.mapin;
offin = p.Results.offin;
sigma = p.Results.sigma;
npoly = p.Results.npoly;
binfac = p.Results.binfac;
dopoly = p.Results.dopoly;
donotch = p.Results.donotch;
filtwidth = p.Results.donotch;
fr1 = p.Results.fr1;

%%
s= size(fr1);
nfr = s(1);
%% make unfiltered map
if domapin == 0
    display('*******making simple unfiltered map*******')
    [dumbmap,dumboff] = fastlinefit_frin(fr1,0,0,1000);
else
    dumbmap = mapin;
    dumboff = offin;
end

if domaskin
    sigmask = maskin;
else
    
    display('*******clip mask*******')
    sigmask = sigmaclipmask(dumbmap,sigma,25);
    im = fitsread_orient('/home/pkorngut/projects/projects2015/Jul/Jul22_2015_starmasks/instmasks/SWIREinstmask.fits');
    sigmask = sigmask.*im;
    
end

quadmap = zeros(1024);
quadmap(1:512,1:512)=1;
quadmap(513:1024,1:512)=2;
quadmap(513:1024,513:1024)=3;
quadmap(1:512,513:1024)=4;


%% assemble the timestreams for the 4 quadrants

for q=1:4
    ts(q).t=0;
end
display('*******assembling the time streams*******')

for i=1:nfr
    display(strcat(num2str(100*i/nfr),'%done',...
        '*******assembling the time streams*******'))

    bb= (squeeze(fr1(i,:,:)) - dumboff).*sigmask;
    for cc = 1:1024
        for rr=1:1024
            for q=1:4
                if quadmap(cc,rr) == q
                    ts(q).t(end+1) =bb(cc,rr);
                end
            end
        end
    end
end
%% some parameters

cp=get_cal_params;
Fs = cp(1).framerate*((1024/2)^2);
dt=1/Fs;
N = length(ts(1).t);
time=(0:(N-1))*dt;
lowtime = 0:dt*binfac:max(time(:));
%% polynomial filtering
if dopoly == 1
    for q=1:4
        display(strcat('binning down quad',num2str(q)))
        sp = find(ts(q).t ~= 0);
        lowy = interp1(time(sp),ts(q).t(sp),lowtime);
        good = find((lowy == lowy));
        ts(q).biny=lowy(good);
        ts(q).binx= lowtime(good);
    end
    
    figure(58)
    for q=1:4
        display(strcat('fitting polynomial quad',num2str(q)))

        p=polyfit(ts(q).binx,ts(q).biny,npoly);
        subplot(2,2,q)
        plot(ts(q).binx,ts(q).biny,'linewidth',2)
        hold on
        plot(ts(q).binx,polyval(p,ts(q).binx),'linewidth',2)
        hold off
        set(gca,'FontSize',22)
        xlabel('Time(s)')
        ylabel('ADU (binned down)')
        v = polyval(p,time);
        ts(q).filt = ts(q).t(:) - v(:);
        v = polyval(p,ts(q).binx);
        ts(q).unfil = ts(q).biny;
        ts(q).biny = ts(q).biny(:) - v(:);
    end
end
%% Look at the power spectra
%totps=0;
%totups=0;
% 
% for q=1:4
%     %subplot(2,2,q)
%     [ts(q).freq,ts(q).ps] = philtimeps(ts(q).binx,ts(q).biny);
%     %loglog(ts(q).freq,(ts(q).ps));
%     totps = totps+ts(q).ps/4;
%     [ts(q).freq,ts(q).ups] = philtimeps(ts(q).binx,ts(q).unfil);
%     totups = totups+ts(q).ups/4;
% 
% end
% 
% figure(112)
% subplot(1,1,1)
% sm=10;
% loglog(ts(1).freq,smooth(totps,sm),'color','red','linewidth',2)
% hold on
% loglog(ts(1).freq,smooth(totups,sm),'color','black','linewidth',2)
% hold off
% set(gca,'FontSize',22)
% xlabel('Frequency (Hz)')
% %k = smooth(totps(:),sm);
% kk = smooth(totups(:),sm);
% 
% axis([.05,100,1e-3,max(kk(:))])
    
%% Notch filtering

if donotch
    Fs = 1./(lowtime(2) - lowtime(1));
    %nus = [9.503,44.82,49.66,63.93,69.56,104,123,163,182,286.6];
    nus = 9.503;
    
    for q=1:4
        dirty = ts(q).biny;
        %display(strcat('Notch filtering quad',num2str(q)))
        
        for tone = 1:numel(nus)
            %if tone*nu1*(1+filtwidth) < Fs/2
                display(strcat('Notch filtering quad',num2str(q),...
                    '--',num2str(nus(tone)),'Hz'))

                d = designfilt('bandstopiir','FilterOrder',2, ...
                    'HalfPowerFrequency1',nus(tone).*(1-filtwidth),...
                    'HalfPowerFrequency2',nus(tone).*(1+filtwidth), ...
                    'DesignMethod','butter','SampleRate',Fs);
                
                clean = filtfilt(d,dirty);
                
                dirty=clean;
            %end
        end
        lowresid = ts(q).biny-clean;
        hiresid = interp1(lowtime,lowresid,time,'spline');

        ts(q).filt = ts(q).filt - hiresid(:);
    end
figure(112)
[f,psn] = philtimeps(lowtime,clean);
hold on
loglog(f,smooth(psn,sm),'color','green','linewidth',2)
hold off
axis([.05,1000,1e-3,max(kk(:))])

end

%%
display('putting filtered timestreams back to frames')

filtfrq=fr1*0;
for q=1:4
count(q).c = 0;
end

for i=1:nfr
    display(strcat(num2str(100*i/nfr),'%done',...
        '*******putting the time streams back*******'))
    
    
    for cc = 1:1024
        for rr=1:1024
            for q=1:4
                if quadmap(cc,rr) == q
                    count(q).c = count(q).c+1;
                    %ts(q).t(end+1) =bb(cc,rr);
                    filtfrq(i,cc,rr) = ts(q).filt(count(q).c);
                end
            end
        end
    end
end
%%
figure(200)
[filtmap,filtoff] = fastlinefit_frin(filtfrq,0,0,1000);
%%
figure(30)
subplot(2,2,1)
imageclip(dumbmap.*sigmask);
title('unfiltered')
subplot(2,2,2)
imageclip(filtmap.*sigmask);
title('filtered')
subplot(2,2,3)
imageclip((dumbmap - filtmap).*sigmask);
title('residual')




return