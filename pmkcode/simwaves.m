np=512;
nfr=20;

cp=get_cal_params;
Fs = cp(1).framerate*((1024/2)^2);
dt=1/Fs;
maxt = dt*(np^2)*nfr;

time = 0:dt:maxt;
nu = 9.503*2;
w = 2*pi*nu;

y = sin(w*time);

%[param]=sine_fit(time,y)
%%
[freq,ps] = philtimeps(time,y);
loglog(freq,ps)
%%

quadmap = zeros(1024);
quadmap(1:512,1:512)=1;
quadmap(513:1024,1:512)=2;
quadmap(513:1024,513:1024)=3;
quadmap(1:512,513:1024)=4;
%%

simfr = zeros(nfr,1024,1024);

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
                    simfr(i,cc,rr) = y(count(q).c);
                end
            end
        end
    end
end
%%

simstripe = fastlinefit_frin(simfr,0);
%%
simstripe =simstripe./max(simstripe(:));

noise = awgn_pk(simstripe*0,1);
%%
load('/home/pkorngut/projects/projects2015/Jul/Jul30_2015_imager_drcombine/fullpipe_Jul30_dr.mat');

f=8;
bigmask = alldat(f).bigmask;
bigmask(10:70,380:520)=0;

%%
nbins=30;

[Cl,l,lth,dCl,binl,dl,Cl2d1] = ...
    get_angular_spec(noise,noise,7.0,nbins,1,ones(1024),'verbose',0,'superbin',0);
lind = NaN;
y1=sqrt(l.*(l+1).*Cl);

[Cl,l,lth,dCl,binl,dl,Cl2d2] = ...
    get_angular_spec(simstripe,simstripe,7.0,nbins,1,ones(1024),'verbose',0,'superbin',0);
lind = NaN;
y2=sqrt(l.*(l+1).*Cl);

[Cl,l,lth,dCl,binl,dl,Cl2d3] = ...
    get_angular_spec(simstripe.*bigmask,simstripe.*bigmask,7.0,nbins,1,ones(1024),'verbose',0,'superbin',0);
lind = NaN;
y3=sqrt(l.*(l+1).*Cl);


[Cl,l,lth,dCl,binl,dl,Cl2d4] = ...
    get_angular_spec(simstripe+noise,simstripe+noise,7.0,nbins,1,ones(1024),'verbose',0,'superbin',0);
lind = NaN;
y4=sqrt(l.*(l+1).*Cl);

%%
subplot(1,2,1)
loglog(l,y2,'linewidth',2)
hold on
loglog(l,y3,'linewidth',2)
%loglog(l,y3)
hold off
set(gca,'FontSize',22,'Xtick',10.^(2:5))
xlabel('l')
ylabel('(l*(l+1)*Cl)^{1/2}')
axis([1e3,1e5,1e-2,10])
legend('Unmasked','Masked')


cax=[0,2e-7];

dx=100;
cx=512;
subplot(2,2,2)
imageclip((Cl2d2));
axis([cx-dx,cx+dx,cx-dx,cx+dx])
caxis(cax)
set(gca,'FontSize',22,'Xtick',10.^(2:5))
title('Unmasked')

subplot(2,2,4)
imageclip((Cl2d3));
axis([cx-dx,cx+dx,cx-dx,cx+dx])
caxis(cax)
set(gca,'FontSize',22,'Xtick',10.^(2:5))
title('Masked')

%%
subplot(1,1,1)
imageclip(simstripe.*bigmask);
