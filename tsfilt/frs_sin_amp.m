function fr_sin=frs_sin_amp(flight,inst,nfr,sin_freq,ampt_poly,phase_in)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%make a time stream with sinwave at sin_freq and amplitude
%follows the polynomial ampt_poly and phase=phase_in
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;
frate=cp(inst).framerate;

Fs = frate*((1024/2)^2);
dt=1/Fs;
N = 512*512*nfr;
time=(0:(N-1))*dt;

ampt_in=polyval(ampt_poly,time);
y =ampt_in.*sin(2*pi*sin_freq.*time+phase_in);

quadmap = zeros(1024);
quadmap(1:512,1:512)=1;
quadmap(513:1024,1:512)=2;
quadmap(513:1024,513:1024)=3;
quadmap(1:512,513:1024)=4;


fr_sin = zeros(nfr,1024,1024);

for q=1:4
count(q).c = 0;
end

for i=1:nfr
    for cc = 1:1024
        for rr=1:1024
            for q=1:4
                if quadmap(cc,rr) == q
                    count(q).c = count(q).c+1;
                    fr_sin(i,cc,rr) = y(count(q).c);
                end
            end
        end
    end
end
return