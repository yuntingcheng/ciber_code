x_arr=1:70;
xlarge_arr=5.5:10:65.5;
xsrc_arr=30.5+10*rand(1,1000);
xsrc_arr=sort(xsrc_arr);
sig=10;
%p_arr=(1./sqrt(2*pi*sig^2)).*exp(-(x_arr-35.5).^2./2./sig.^2);
%%
figure
setwinsize(gcf,2000,200)

pstack_arr=zeros(size(xlarge_arr));
stamp=zeros(1,61);
count=0;
for isrc=1:numel(xsrc_arr)
    xsrc=xsrc_arr(isrc);
    p_arr=(1./sqrt(2*pi*sig^2)).*exp(-(x_arr-xsrc).^2./2./sig.^2);
    
    if rem(isrc,300)==1
        count=count+1;
        subplot(1,4,1)
        plot(x_arr,p_arr,'o-','linewidth',2,'color',get_color(count));hold on
        title('input PSF')
    end
    
    plarge_arr=zeros(size(xlarge_arr));
    psmall_arr=zeros(size(x_arr));
    for ipix=1:numel(xlarge_arr)
        plarge_arr(ipix)=mean(p_arr(ipix*10-9:ipix*10));
        psmall_arr(ipix*10-9:ipix*10)=plarge_arr(ipix);
    end
    
    if rem(isrc,300)==1
        subplot(1,4,2)
        plot(xlarge_arr,plarge_arr,'o-','linewidth',2,...
            'color',get_color(count));hold on
        title('bin to large pixel')
    end
    pstack_arr=pstack_arr+plarge_arr./numel(xsrc_arr);
    
    stamp=stamp+psmall_arr(round(xsrc-30:xsrc+30))./numel(xsrc_arr);
    if rem(isrc,300)==1
        subplot(1,4,3)
        plot(-30:30,psmall_arr(round(xsrc-30:xsrc+30)),...
            'linewidth',2,'color',get_color(count));hold on
        title('small pixel stamps')
    end
    
end

subplot(1,4,4)
plot(-30:30,stamp,'linewidth',2);hold on
plot(-30:30,(1./sqrt(2*pi*sig^2)).*exp(-((1:61)-31).^2./2./sig.^2),'linewidth',2);
legend({'stack PSF','true PSF'})
title('stack 1000 random sources')

savedir=(strcat('/Users/ytcheng/ciber/doc/20171130_psfstack/plots/'));

savename=strcat(savedir,'stack_test1D');
print(savename,'-dpng');%close