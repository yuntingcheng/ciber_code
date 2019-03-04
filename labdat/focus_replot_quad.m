%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For some TM2 data, the top right quadrant are dead,
%which makes the map plotted in slopemap has a high color range and 
%hard to see the structure. So this code generate plot map again 
%by offset the top right quadrant to the median of rest of the map.
%The results is stored kin the same folder with name ...._fix.png.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
date='03-13-2013';
band=2;
savedatdir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/',...
            'Focus/slopedata/',date,'/TM',num2str(band),'/'); 
saveplotdir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/',...
            'Focus/slopemap/',date,'/TM',num2str(band),'/'); 

scantime=dir(strcat(savedatdir,'*framedat.mat'));
for i=1:numel(scantime)
time=scantime(i).name(1:8);
load(strcat(savedatdir,scantime(i).name),'framedat');
map=squeeze(framedat.map_arr(end,:,:));
goodpix1=map(1:512,1:512);goodpix1=goodpix1(:);
goodpix2=map(1:512,513:1024);goodpix2=goodpix2(:);
goodpix3=map(513:1024,1:512);goodpix3=goodpix3(:);
goodpix=[goodpix1 goodpix2 goodpix3];goodpix=goodpix(:);
map(513:1024,513:1024)=median(goodpix);

imageclip(map);
plotname=strcat(saveplotdir,'C',num2str(band),'_',date,'_',time,'_fix');
print(plotname,'-dpng');close
end