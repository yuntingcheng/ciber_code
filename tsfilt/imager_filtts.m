function [filtfrq] = imager_filtts(fr1,offin,maskin,varargin)
 
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('fr1');
p.addRequired('offin',@isnumeric); 
p.addRequired('maskin',@isnumeric);
p.addOptional('flight',40030,@isnumeric);
p.addOptional('inst',1,@isnumeric);
p.addOptional('npoly',3,@isnumeric); 
p.addOptional('binfac',100,@isnumeric); 
p.addOptional('dopoly',0,@isnumeric); 
p.addOptional('donotch',0,@isnumeric); 
p.addOptional('dosinfilt',1,@isnumeric);
p.addOptional('dorangesinfilt',0,@isnumeric); 
p.addOptional('doampvtime',1,@isnumeric); 
p.addOptional('np_ampvtime',8,@isnumeric); 
p.addOptional('sin_freq',9.503,@isnumeric);
p.addOptional('sin_npoly',4,@isnumeric); 
p.addOptional('notch_nus',9.503,@isnumeric);
p.addOptional('notch_width',0.5,@isnumeric); 
p.addOptional('notch_iter',1,@isnumeric); 
p.addOptional('verbose',0,@isnumeric); 
p.addOptional('makeplot',0,@isnumeric); 
p.addOptional('saveplot',0,@isnumeric);
p.addOptional('plotdir','',@isstr);
p.parse(fr1,offin,maskin,varargin{:});

fr1 = p.Results.fr1;
offin = p.Results.offin;
maskin = p.Results.maskin;
flight = p.Results.flight;
inst = p.Results.inst;
npoly = p.Results.npoly;
binfac = p.Results.binfac;
dopoly = p.Results.dopoly;
donotch = p.Results.donotch;
dosinfilt = p.Results.dosinfilt;
dorangesinfilt = p.Results.dorangesinfilt;
doampvtime = p.Results.doampvtime;
np_ampvtime = p.Results.np_ampvtime;
sin_freq = p.Results.sin_freq;
sin_npoly = p.Results.sin_npoly;
notch_nus = p.Results.notch_nus;
notch_width = p.Results.notch_width;
notch_iter = p.Results.notch_iter;
verbose = p.Results.verbose;
makeplot = p.Results.makeplot;
saveplot = p.Results.saveplot;
plotdir = p.Results.plotdir;

%%
s= size(fr1);
nfr = s(1);
%% assemble the timestreams for the 4 quadrants
dumboff = offin;
sigmask = maskin;

quadmap = zeros(1024);
quadmap(1:512,1:512)=1;
quadmap(513:1024,1:512)=2;
quadmap(513:1024,513:1024)=3;
quadmap(1:512,513:1024)=4;

for q=1:4
    ts(q).t=[];
end

for i=1:nfr
    
    if verbose
        display(strcat(num2str(100*i/nfr),'%done',...
            '*******assembling the time streams*******'))
    end
    
    bb= (squeeze(fr1(i,:,:)) - dumboff).*sigmask;
    
    %%%% this is already the correct scan direction, don't need to flip %%%
    %bb(513:1024,1:512)=flip(bb(513:1024,1:512),1);
    %bb(1:512,513:1024)=flip(bb(1:512,513:1024),2);
    %bb(513:1024,513:1024)=flip(bb(513:1024,513:1024),1);
    %bb(513:1024,513:1024)=flip(bb(513:1024,513:1024),2);
    
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

for q=1:4
    ts(q).t=ts(q).t(:);
end
%% some parameters
cp=get_cal_params('flight',flight);
Fs = cp(inst).framerate*((1024/2)^2);
dt=1/Fs;
N = length(ts(1).t);
time=(0:(N-1))*dt;
lowtime = 0:dt*binfac:max(time(:));
%% plot the time stream without filtering
if makeplot
    figure
    setwinsize(gcf,800,600)
    for q=1:4
     [freq_arr,ups_arr] = philtimeps(time,ts(q).t);
     subplot(2,2,q)
     loglog(freq_arr,ups_arr);hold on
     xlabel('Frequency (Hz)')
     drawnow
    end
end
%% polynomial filtering
if dopoly == 1
    for q=1:4
        sp = find(ts(q).t ~= 0);
        lowy = interp1(time(sp),ts(q).t(sp),lowtime);
        good = find((lowy == lowy));
        ts(q).biny=lowy(good);
        ts(q).binx= lowtime(good);
        ts(q).good=good;
    end
    
    if makeplot
        figure
        setwinsize(gcf,1000,1000)
    end
    
    for q=1:4
        %display(strcat('fitting polynomial quad',num2str(q)))

        p=polyfit(ts(q).binx,ts(q).biny,npoly);
        if makeplot
            subplot(2,2,q)
            plot(ts(q).binx,ts(q).biny,'linewidth',2)
            hold on
            plot(ts(q).binx,polyval(p,ts(q).binx),'linewidth',2)
            hold off
            set(gca,'FontSize',22)
            xlabel('Time(s)')
            ylabel('ADU (binned down)')
            drawnow
        end
        v = polyval(p,time);
        ts(q).filt = ts(q).t(:) - v(:);
        v = polyval(p,ts(q).binx);
        ts(q).unfil = ts(q).biny;
        ts(q).biny = ts(q).biny(:) - v(:);
    end

%%% Look at the power spectra

if makeplot
    figure
    setwinsize(gcf,1000,1000)
end

 for q=1:4
     %subplot(2,2,q)
     [ts(q).freq,ts(q).ps] = philtimeps(ts(q).binx,ts(q).biny);
     %loglog(ts(q).freq,(ts(q).ps));
     %totps = totps+ts(q).ps/4;
     totps=ts(q).ps;
     [ts(q).freq,ts(q).ups] = philtimeps(ts(q).binx,ts(q).unfil);
     %totups = totups+ts(q).ups/4; 
     totups=ts(q).ups;
if makeplot
     subplot(2,2,q)
     sm=10;
     loglog(ts(q).freq,smooth(totps,sm),'color','red','linewidth',1)
     hold on
     loglog(ts(q).freq,smooth(totups,sm),'color','black','linewidth',1)
     hold off
     %set(gca,'FontSize',22)
     xlabel('Frequency (Hz)')
     %k = smooth(totps(:),sm);
     kk = smooth(totups(:),sm);
     axis([.05,100,1e-3,max(kk(:))])
     drawnow
end
end

 
else
    for q=1:4
    ts(q).filt = ts(q).t(:);
    end
end
%% Notch filtering
if donotch
    npad1=round(numel(ts(q).t)*0.2);
    npad2=round(numel(ts(q).t)*0.1);
    %Fs = 1./(lowtime(2) - lowtime(1));
    Fs = 1./dt;
    %notch_nus = [9.503,44.82,49.66,63.93,69.56,104,123,163,182,286.6];
    %notch_nus=[9.503,286,286*2,286*3,286*4];
    
    for q=1:4
        
        dirty=[zeros(1,npad1) ts(q).t' zeros(1,npad2)];
        dirty=dirty(:);
        
        for tone = 1:numel(notch_nus)
            %if tone*nu1*(1+filtwidth) < Fs/2
            if verbose
                display(strcat('Notch filtering quad',num2str(q),...
                    '--',num2str(notch_nus(tone)),'Hz'))
            end
            %d = designfilt('bandstopiir','FilterOrder',2, ...
            %    'HalfPowerFrequency1',notch_nus(tone).*(1-filtwidth),...
            %    'HalfPowerFrequency2',notch_nus(tone).*(1+filtwidth), ...
            %    'DesignMethod','butter','SampleRate',Fs);


            for inotch_iter=1:notch_iter
            d = designfilt('bandstopiir','FilterOrder',2, ...
                'HalfPowerFrequency1',notch_nus(tone)-notch_width(tone),...
                'HalfPowerFrequency2',notch_nus(tone)+notch_width(tone), ...
                'DesignMethod','butter','SampleRate',Fs);
            clean = filtfilt(d,dirty);  
            dirty=clean;
            end
            
        end
        %lowresid = ts(q).biny-clean;
        %hiresid = interp1(lowtime(good),lowresid(:),time,'spline');
        %ts(q).filt = ts(q).filt - hiresid(:);
        
        ts(q).filt=clean(npad1+1:end-npad2);
    end


if makeplot
    for tone = 1:numel(notch_nus)
        figure
        setwinsize(gcf,800,600)
        for q=1:4
         [freq_arr,ups_arr] = philtimeps(time,ts(q).t);
         [~,ps_arr] = philtimeps(time,ts(q).filt);
         subplot(2,2,q)
         loglog(freq_arr,ups_arr,'k');hold on
         loglog(freq_arr,ps_arr,'r');hold off
         xlabel('Frequency (Hz)')
         xlim([notch_nus(tone)-10*notch_width(tone),...
             notch_nus(tone)+10*notch_width(tone)])
         drawnow
        end
    end
end
%%
end
%% Sinewave filtering
%{
for q=1:4
ts(q).filt = ts(q).t(:);
end
%}
if dosinfilt
    if verbose
        display('Filtering out Sin Wave')
    end
    
    for q=1:4
        [filtts,stack,stackfit,bigsinparams,tchu] = ...
            pulsar_filtsin(time(:),ts(q).filt(:),sin_freq);
        if doampvtime == 0
            ts(q).filt = filtts(:);
        else
            np=8;
            top = round(numel(time)/np);
            allts = (1:np)*0;
            allamps = allts;
            for p = 1:np
                tp = time((p-1)*top +1:p*top);
                mp = ts(q).t((p-1)*top +1:p*top);
                [~,~,~,sinparams] = pulsar_filtsin(tp(:),mp(:),sin_freq);
                allts(p) = mean(tp(:));
                allamps(p) = sinparams.amp;
            end
            po=polyfit(allts,allamps,sin_npoly);
            ampt=polyval(po,time);
            bigmodel = ampt.*sin(2*pi*sin_freq.*time+bigsinparams.phase);
            ts(q).filt = ts(q).filt(:)-bigmodel(:);
            sinampdat(q).allts=allts;
            sinampdat(q).allamps=allamps;
            sinampdat(q).ampt=ampt;
            sinampdat(q).bigamp=bigsinparams.amp;
            sinampdat(q).bigmodel=bigmodel; 
            
            sinpar(q).po=po;
            sinpar(q).phase=bigsinparams.phase;
        end
    end
%%   
if makeplot 
    %%
    figure
    setwinsize(gcf,800,600)
    for q=1:4
    subplot(2,2,q)
    plot(sinampdat(q).allts,sinampdat(q).allamps,'bo');hold on
    plot(time,sinampdat(q).ampt,'b');
    %hline = refline([0 sinampdat(q).bigamp]);
    %hline.Color = 'k';
    xlabel('t(s)');
    legend({'data','polyfit'},'location','northwest');
    drawnow
    end
    if saveplot
    savename=strcat(plotdir,'sinampfit');
    print(savename,'-dpng');close
    end
    %%
     figure
     setwinsize(gcf,800,400)
     [freq_arr,ups_arr] = philtimeps(time,ts(q).t);
     [~,ps_arr] = philtimeps(time,ts(q).filt);
     subplot(1,2,1)
     loglog(freq_arr,ups_arr);hold on
     loglog(freq_arr,ps_arr);
     xlabel('Frequency (Hz)')
     xlim([min(freq_arr),max(freq_arr)])
     subplot(1,2,2)
     loglog(freq_arr,ups_arr);hold on
     loglog(freq_arr,ps_arr);
     xlabel('Frequency (Hz)')
     legend({'unfilt','filt'})
     xlim([8.5,10.5])
     drawnow
     %%
    if saveplot
    savename=strcat(plotdir,'sinfiltamp_ps');
    print(savename,'-dpng');close
    end
    
     figure
     setwinsize(gcf,800,400)
     subplot(1,2,1)
     plot(time,ts(q).t);hold on
     plot(time,sinampdat(q).bigmodel);
     xlabel('t (s)')
     subplot(1,2,2)
     plot(time,ts(q).t);hold on
     plot(time,sinampdat(q).bigmodel);     
     ylim([min(sinampdat(q).bigmodel)*10,max(sinampdat(q).bigmodel)*10])
     xlim([time(1),time(end)])
     xlabel('t (s)')
     legend({'unfilt time stream','filt sin wave'})
     drawnow
    if saveplot
    savename=strcat(plotdir,'sinfiltamp_ts');
    print(savename,'-dpng');close
    end   
    %%
end


end
%% Range of Sinewave filtering
np=np_ampvtime;
top = round(numel(time)/np);

[freq_arr,ps_arr]=philtimeps(time,ts(1).t);

if dorangesinfilt
freqs=freq_arr(find(freq_arr<9.6 & freq_arr>9.4));
tic
for q=1:4
dirty = ts(q).filt;
for iter=1:2
for f=1:numel(freqs)
    if verbose
    disp(sprintf('quad%d--iter%d--freq%d/%d(%.3fHz)--%.2f min',...
                    q,iter,f,numel(freqs),freqs(f),toc/60));
    end

    [clean,~,~,bigsinparams] = pulsar_filtsin(time,dirty,freqs(f));
    if doampvtime==0
        dirty=clean;
    else
        allts = (1:np)*0;
        allamps = allts;
        for p = 1:np
            tp = time((p-1)*top +1:p*top);
            mp = dirty((p-1)*top +1:p*top);
            [~,~,~,sinparams] = pulsar_filtsin(tp(:),mp(:),freqs(f));
            allts(p) = mean(tp(:));
            allamps(p) = sinparams.amp;
        end
        po=polyfit(allts,allamps,5);
        ampt=polyval(po,time);
        bigmodel = ampt.*sin(2*pi*freqs(f).*time+bigsinparams.phase);
        clean = dirty-bigmodel(:);
        dirty=clean;
    end
end
end
ts(q).filt = clean;
end


if makeplot
    figure
    setwinsize(gcf,800,400)
    for q=1%1:4
     [~,ups_arr] = philtimeps(time,ts(q).t);
     [~,ps_arr] = philtimeps(time,ts(q).filt);
     subplot(1,2,1)
     loglog(freq_arr,ups_arr);hold on
     loglog(freq_arr,ps_arr);hold off
     xlabel('Frequency (Hz)')
     xlim([min(freq_arr),max(freq_arr)])
     subplot(1,2,2)
     loglog(freq_arr,ups_arr);hold on
     loglog(freq_arr,ps_arr);hold off
     xlabel('Frequency (Hz)')
     legend({'unfilt','filt'})
     xlim([9.2,9.8])
     drawnow
    end
    if saveplot
    savename=strcat(plotdir,'sinfiltrange');
    print(savename,'-dpng');close
    end
end


end
%% putting the time streams back to frames

filtfrq=fr1*0;
for q=1:4
count(q).c = 0;
end

for i=1:nfr
    
    if verbose
    display(strcat(num2str(100*i/nfr),'%done',...
        '*******putting the time streams back*******'))
    end
    
    
    for cc = 1:1024
        for rr=1:1024
            for q=1:4
                if quadmap(cc,rr) == q
                    count(q).c = count(q).c+1;
                    %ts(q).t(end+1) =bb(cc,rr);
                    if maskin(cc,rr)
                        filtfrq(i,cc,rr) = ts(q).filt(count(q).c);
                    else
                        filtfrq(i,cc,rr) = squeeze(fr1(i,cc,rr,:));
                    end
                end
            end
        end
    end
end

return