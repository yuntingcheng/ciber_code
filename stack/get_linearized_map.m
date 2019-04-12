function get_linearized_map(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get linearized mask for stacking
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

loaddir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
load(strcat(loaddir,'FFdat'),'FFdat');
load(strcat(loaddir,'maskdat'),'maskdat');

mask_inst = maskdat.mask(4).mask_inst;

for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    [fr] = get_data_frames(inst,dt.name,'flight',flight,'verbose',0);
    fr=fr(3:end,:,:);
    
    %%% mask pixel having crazy time stream
    dfr=fr(2:end,:,:)-fr(1:end-1,:,:);

    fixfr = fr;
    shitmask = ones(1024);
    shitcut=5e4;
    count1=0;count2=0;count3=0;
    for cc=1:1024
    for rr=1:1024

        d1 = dfr(:,cc,rr);
        sp = find(d1 > shitcut);

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
        %plot(fixfr(:,cc,rr) - fixfr(1,cc,rr));%%%%
    end
    end
    mask=mask_inst.*shitmask;
    
    %%% line fit
    [long,off] = linfit_map(fixfr,'verbose',0);
    [short] = linfit_map(fixfr(1:4,:,:),'verbose',0);

    longcal = long./FFdat(ifield).FFunholy;
    shortcal = short./FFdat(ifield).FFunholy;

    %%% replace saturated pixel with short line fit
    qq= squeeze(fixfr(end,:,:));
    qqs = qq -off;
    sp = find(qqs < -5000);
    fixmap = longcal;
    fixmap(sp) = shortcal(sp);
    
    %%% for elat30, make another map only last 10 frames
    % since elat30 field integration is very unstable in ~ first half
    if ifield == 5
        fixfr = fixfr(end-9:end,:,:);
        [long1,off1] = linfit_map(fixfr,'verbose',0);
        [short1] = linfit_map(fixfr(1:4,:,:),'verbose',0);

        longcal1 = long1./FFdat(ifield).FFunholy;
        shortcal1 = short1./FFdat(ifield).FFunholy;

        qq= squeeze(fixfr(end,:,:));
        qqs = qq -off1;
        sp = find(qqs < -2000);
        fixmap1 = longcal1;
        fixmap1(sp) = shortcal1(sp);
    end
    
    %%% mask negative flux pixel (positive ADU/fr)
    negmask=ones(1024);
    negmask(find(fixmap>0))=0;
    mask = mask.*negmask;

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
    
    stackmapdat(ifield).rawmap = long;
    stackmapdat(ifield).map = fixmap;
    stackmapdat(ifield).mask = mask;
    if ifield == 5
        stackmapdat(ifield).rawmap_last10 = long1;
        stackmapdat(ifield).map_last10 = fixmap1;
    end
    
end

savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
save(sprintf('%s/stackmapdat',savedir),'stackmapdat');

return