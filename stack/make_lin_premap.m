flight=40030;
inst=2;
mypaths=get_paths(flight);
cp=get_cal_params('flight',flight);
cal=cp(inst).apf2eps.*cp(inst).eps2nWpm2ps;

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'FFdat'),'FFdat');
load(strcat(loaddir,'maskdat'),'maskdat');

for ifield=4:8

    %%% get the frames
    dt=get_dark_times(flight,inst,ifield);
    [fr_all] = get_data_frames(inst,dt.name,'flight',flight,'verbose',0);

    fr=fr_all(3:end,:,:);

    %%% mask pixel having crazy time stream
    dfr=fr(2:end,:,:)-fr(1:end-1,:,:);

    fixfr = fr;
    shitmask = ones(1024);
    shitcut1=1e4;
    shitcut2=5e4;
    count1=0;count2=0;count3=0;
    for cc=1:1024
    for rr=1:1024

        d1 = dfr(:,cc,rr);
        sp = find(d1 > shitcut2);

        if numel(sp)>1
            shitmask(cc,rr)=0;
            count1=count1+1;
        elseif numel(sp)==1
            f1 = fr(:,cc,rr);
            f1(sp+1:end) = f1(sp+1:end) - 2^16;
            fixfr(:,cc,rr)=f1;
            count2=count2+1;
        elseif numel(find(abs(d1)>median(abs(d1))*100))>0
            shitmask(cc,rr)=0;
            count3=count3+1;
        end
    end
    end
    mask=maskdat.mask(ifield).mask_inst.*shitmask;

    %%% line fit and FF/cal correction

    [long,off] = linfit_map(fixfr,'verbose',0);
    [short] = linfit_map(fixfr(1:4,:,:),'verbose',0);

    FF=FFdat(ifield).FF;
    FFmask=ones(1024);FFmask(find(FF==0))=0;
    flat = unholy_map(FF,FFmask,200);
    longcal = cal*long./flat;
    shortcal = cal*short./flat;

    %%% linearize
    lastfr= squeeze(fixfr(end,:,:));
    lastfroff =lastfr -off;
    sp = find(lastfroff < -5000);

    map = longcal;
    map(sp) = shortcal(sp);

    mask(find(map<0))=0;

    %%% hand mask bad pixels
    if inst==1 & ifield==8
    mask=circular_mask(442,40,40,mask);
    mask=circular_mask(742,342,10,mask);
    end

    if inst==2 & ifield==8
    mask=circular_mask(450,395,10,mask);
    mask=circular_mask(0,800,60,mask);
    end

    if inst==2 & ifield==5
    mask=elliptical_mask(215,220,60,15,80,mask);
    mask=elliptical_mask(710,90,80,20,60,mask);
    end
    
    %%% save the map and mask
    savedir=(strcat('/Users/ytcheng/ciber/doc/20171018_stackihl/stackmaps/TM',...
        num2str(inst),'/'));

    fits_write(strcat(savedir,'maps/',dt.name,'_map'),map);
    fits_write(strcat(savedir,'masks/',dt.name,'_mask'),mask);

    fits_write(strcat(savedir,'maps/',dt.name,'_mapA'),map(1:512,1:512));
    fits_write(strcat(savedir,'maps/',dt.name,'_mapB'),map(513:1024,1:512));
    fits_write(strcat(savedir,'maps/',dt.name,'_mapC'),map(1:512,513:1024));
    fits_write(strcat(savedir,'maps/',dt.name,'_mapD'),map(513:1024,513:1024));

end
