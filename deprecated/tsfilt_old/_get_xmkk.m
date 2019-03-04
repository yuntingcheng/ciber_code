savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
load(strcat(savedirfull,'fwfulldat'),'fwfulldat');
pixscale=7;

count=0;
[~,~,~,~,binl] = get_angular_spec(randn(1024),randn(1024),pixscale);
for ifield=[8,7,6,5,4,2,1]
    for jfield=[8,7,6,5,4,2,1]
        if ifield>jfield
            count=count+1;
            disp(sprintf('count=%d,[i,j]=[%d,%d]',count,ifield,jfield))
            mask=fwfulldat(ifield).mask.*fwfulldat(jfield).mask;
            std_noisei=fwfulldat(ifield).std_fCl2d;
            weighti=(fftshift(fftshift(1./std_noisei)))';
            std_noisej=fwfulldat(jfield).std_fCl2d;
            weightj=(fftshift(fftshift(1./std_noisej)))';
            weight=(weighti+weightj)./2;
            Mkk =get_mkk_sim(mask,pixscale,binl,100,...
                numel(binl),1,weight,0,NaN);
            xmkk(count).field=[ifield jfield];
            xmkk(count).Mkk=Mkk;
        end
    end
end

savedir='/Users/ytcheng/ciber/doc/20170325_alldat/TM1/';
save(strcat(savedir,'xmkk'),'xmkk');