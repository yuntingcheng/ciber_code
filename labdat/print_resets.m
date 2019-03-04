%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For the given band, scan trough the dark data and find resets.
% print out the reset frames number which 90% up of pixel find 
% it as a reset.
%
% Output (print out):
% - time: dark data time label ex: 13-53-13
% - reset: reset frames number
% - fraction: fraction of pix find this frame is a reset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
band=2;
flight=36265;

switch flight
case {36265,36277}
darkdir=strcat('/Users/ytcheng/ciber/data/Dark_Raw/362xx/TM',...
    num2str(band),'/');
case 40030
darkdir=strcat('/Users/ytcheng/ciber/data/Dark_Raw/40030/TM',...
    num2str(band),'/');
end    
%% find resets and print out 
scanfile=dir(strcat(darkdir,'*0001.fits'));
for i=1:numel(scanfile)
    datetime=scanfile(i).name(8:26);
    files = dir(strcat(darkdir,'*',datetime,'*')); 
    
    %%%%%%%%read in frames data%%%%%%%%%%
    frames=zeros(numel(files),1024,1024);
    for ifr=1:numel(files)
    infile = strcat(darkdir,files(ifr).name);
    frame = imrotate(fitsread(infile),270);    
    frames(ifr,:,:)=frame; 
    end
    
    %%%%%%%%get resets stats%%%%%%%%%%%%%
    reset_arr=zeros(numel(files),1024,1024);
    for ir=1:1024
    for ic=1:1024
    ts=frames(:,ir,ic);
    reset_arr(:,ir,ic)=get_resets(ts);
    end
    end
    rstat_arr=squeeze(sum(reset_arr,3));
    rstat_arr=squeeze(sum(rstat_arr,2));
    rstat_arr=rstat_arr./1024./1024;

    %%%%%%%print results%%%%%%%%%%%%%%%%%
    pr=sprintf('%s, Nfr=%d',datetime,numel(files));disp(pr);
    find(rstat_arr>0.9)
    rstat_arr(find(rstat_arr>0.9))
end
