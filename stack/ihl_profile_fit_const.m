%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot ihl stacking maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
flight=40030;
inst=1;
ifield=8;
npix=800;

dt=get_dark_times(flight,inst,ifield);

savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/ihl_const/TM',...
    num2str(inst),'/'));
quad_arr=['A','B','C','D'];
%% get the stacking map
m_arr=13:17;
for im=1:numel(m_arr)
m=m_arr(im);
%%% read in 2D maps
for itype=1:2
    if itype==1
        type='ciber_2m';
    else
        type='tmass_2m';
    end
      loaddir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
    num2str(inst),'/',char(type),'/'));
  
    stackgtot=zeros(npix*2+1);
    mstackgtot=zeros(npix*2+1);
    for iquad=1:4
        quad=quad_arr(iquad);
        if itype==1
            stackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
                '_stamperg',num2str(m),'.fits'));
        else
            stackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
                '_stamperg',num2str(m),'.fits'));            
        end
        
        mstackg=fitsread(strcat(loaddir,dt.name,'_',quad,...
            '_maskstamperg',num2str(m),'.fits'));

        stackgtot=stackgtot+stackg;
        mstackgtot=mstackgtot+mstackg;  
    end
    if itype==1
        stackgcb=stackgtot./mstackgtot;
        stackgcb=stackgcb./stackgcb(npix+1,npix+1);
    else
        stackg2m=stackgtot./mstackgtot;
        stackg2m=stackg2m./stackg2m(npix+1,npix+1);
    end  
    
    
end
diffmap=stackgcb-stackg2m;

%%%%% take the 1D profile %%%%%
profile = radial_prof(diffmap,ones(2*npix+1),npix+1,npix+1,1,200);
r_arr=profile.r*0.7;
sp=find(r_arr>60);
prof_arr=profile.prof;
err_arr=profile.err;
off=sum(prof_arr(sp)./err_arr(sp).^2)/sum(1./err_arr(sp).^2);

figure
xlimit=[10,800];
ylimit=[off-3*std(prof_arr(sp)),off+10*std(prof_arr(sp))];
errorbar(r_arr,prof_arr,err_arr,'b','DisplayName','CIBER-2MASS');hold on
plot(xlimit,[off,off],'r','DisplayName',sprintf('offset=%.2e',off));
plot([60,60],ylimit,'k','DisplayName','fit offset (r>1 arcmin)');
title(strcat(num2str(m),'<m<',num2str(m+1)),'FontSize', 20)
xlabel('arcsec','FontSize', 20)
ylabel('CIBER stack - 2MASS stack','FontSize', 20)

h=legend('show','Location','northeast');
set(h,'fontsize',20)
legend boxoff

set(gca,'XScale','log')
xlim(xlimit)
ylim(ylimit)

savename=strcat(savedir,dt.name,'_',num2str(m),'_diffprof');
print(savename,'-dpng');%close

offsetdat(im).m=m;
offsetdat(im).off=off;
end

save(strcat(savedir,'offsetdat'),'offsetdat');
