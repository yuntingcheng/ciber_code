band=2;

datadir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/fitTM'...
    ,num2str(band),'/'); 

scanfile=dir(sprintf('%snfr2*',datadir));
disp(sprintf('TM%d,nfr2 has %d files',band,numel(scanfile)));
%%
for ifile=10%1:numel(scanfile)

dataname=sprintf('%s%s',datadir,scanfile(ifile).name);
load(dataname,'fitdat');

for imap=1:numel(fitdat.avg_arr)
diffmap=squeeze(fitdat.diffmap_arr(imap,:,:));
figure
imageclip(diffmap);
title(sprintf('ifile=%d,var=%.4f',ifile,fitdat.var_arr(imap)));
disp(sprintf('imap=%d,std=%.4f',imap,...
            std(diffmap(find(diffmap)))));
end
end
