flight=40030;
inst=1;
ifield=8;
mypaths=get_paths(flight);
dt=get_dark_times(flight,inst,ifield);
savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/'));
m_arr=10:1:17;

%% masking parameters (mask is saved, don't run it everytime)
%{
alpha=-9.7*1.2;
beta=159*1.2;

for m_max=m_arr
    
pscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','PSC');
xscmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSC');
xscrmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,'catname','XSCrej');
pscrmask=make_strmask_2m(flight,inst,ifield,alpha,beta,m_max,...
'catname','PSCrej','rel','A');
strmask=pscmask.*xscmask.*xscrmask.*pscrmask;

fits_write(strcat(savedir,'masks2m/',dt.name,'_strmask',num2str(m_max)),strmask);
disp(sprintf('strmask within m=%d save',m_max));
end
%}
%% get inst mask and map
mask=fits_read(strcat(savedir,'masks/',dt.name,'_mask.fits'));
map=fits_read(strcat(savedir,'maps/',dt.name,'_map.fits'));

% correct for the cal factor error
map=map./0.383;

srcmapdir=strcat('/Users/ytcheng/ciber/doc/20170617_Stacking/srcmaps/TM',...
    num2str(inst),'/');
quad_arr=['A','B','C','D'];
tmmap=zeros(1024);
for iquad=1:4
    quad=quad_arr(iquad);
    stmmaps=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmaps_2m.fits'));
    stmmapg=fits_read(strcat(srcmapdir,dt.name,'_',quad,'_srcmapg_2m.fits'));
    stmmap=stmmaps+stmmapg;
    
    if iquad==1
        tmmap(1:512,1:512)=stmmap;
    elseif iquad==2
        tmmap(513:1024,1:512)=stmmap;
    elseif iquad==3
        tmmap(1:512,513:1024)=stmmap;
    else
        tmmap(513:1024,513:1024)=stmmap;
    end
end

%% get diff mask mask to 10 mag
maskin=mask;
clipmin=0;
clipmax=1.2298e+05;% from make_magstrmask.m UK process
%%

for m_max=m_arr
%%% inst mask * strmask
strmask=fits_read(strcat(savedir,'masks2m/',dt.name,'_strmask',...
    num2str(m_max),'.fits'));   

totmask=maskin.*strmask;
sp=find(map>clipmax | map<clipmin);
totmask(sp)=0;

meancb=median(map(find(totmask)));
mean2m=median(tmmap(find(totmask)));

mmapcb=(map-meancb+mean2m).*totmask;
mmap2m=tmmap.*totmask;

fits_write(strcat(savedir,'maps2m/',dt.name,'_cbmap',num2str(m_max),'A'),...
    mmapcb(1:512,1:512));
fits_write(strcat(savedir,'maps2m/',dt.name,'_cbmap',num2str(m_max),'B'),...
    mmapcb(513:1024,1:512));
fits_write(strcat(savedir,'maps2m/',dt.name,'_cbmap',num2str(m_max),'C'),...
    mmapcb(1:512,513:1024));
fits_write(strcat(savedir,'maps2m/',dt.name,'_cbmap',num2str(m_max),'D'),...
    mmapcb(513:1024,513:1024));

fits_write(strcat(savedir,'maps2m/',dt.name,'_2mmap',num2str(m_max),'A'),...
    mmap2m(1:512,1:512));
fits_write(strcat(savedir,'maps2m/',dt.name,'_2mmap',num2str(m_max),'B'),...
    mmap2m(513:1024,1:512));
fits_write(strcat(savedir,'maps2m/',dt.name,'_2mmap',num2str(m_max),'C'),...
    mmap2m(1:512,513:1024));
fits_write(strcat(savedir,'maps2m/',dt.name,'_2mmap',num2str(m_max),'D'),...
    mmap2m(513:1024,513:1024));

fits_write(strcat(savedir,'masks2m/',dt.name,'_totmask',num2str(m_max),'A'),...
    totmask(1:512,1:512));
fits_write(strcat(savedir,'masks2m/',dt.name,'_totmask',num2str(m_max),'B'),...
    totmask(513:1024,1:512));
fits_write(strcat(savedir,'masks2m/',dt.name,'_totmask',num2str(m_max),'C'),...
    totmask(1:512,513:1024));
fits_write(strcat(savedir,'masks2m/',dt.name,'_totmask',num2str(m_max),'D'),...
    totmask(513:1024,513:1024));

end


