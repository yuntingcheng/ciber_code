function get_bl_old(flight,inst)
mypaths=get_paths(flight);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');

nsim=50;
pixsize=0.7;

for ifield=4:8
    dt=get_dark_times(flight,inst,ifield);
    name=dt.name;
    if inst==1
        foldname='/Users/ytcheng/ciber/data/iband_psf_quads/';
    end

    a=fitsread(strcat(foldname,name,'_A.fits'));
    b=fitsread(strcat(foldname,name,'_B.fits'));
    c=fitsread(strcat(foldname,name,'_C.fits'));
    d=fitsread(strcat(foldname,name,'_D.fits'));

    ah=fitsread(strcat(foldname,name,'_Ahitmap.fits'));
    bh=fitsread(strcat(foldname,name,'_Bhitmap.fits'));
    ch=fitsread(strcat(foldname,name,'_Chitmap.fits'));
    dh=fitsread(strcat(foldname,name,'_Dhitmap.fits'));

    psfmap=(a.*ah+b.*bh+c.*ch+d.*dh)./(ah+bh+ch+dh);
    rad = make_radius_map(psfmap,401,401);
    psfmap=psfmap-mean(psfmap(find(rad>100)));
    psfmap=psfmap./sum(psfmap(:));
    
    bigpsf=get_coarse_map(psfmap,10);
    
    bldat.psfmap(ifield).psfmap=psfmap;
    bldat.psfmap(ifield).bigpsf=bigpsf;
    
    bl_arr=zeros(nsim,29);
    for isim=1:nsim
        disp(sprintf('ifield%d,isim%d',ifield,isim));
        map=randn(1024);map=map-mean(map(:));
        C=conv2(map,psfmap(301:501,301:501));
        % use the smaller psf to speed up the code,
        % I tested the results is consistent with using full psf
        Nout=(size(C,1)-size(map,1))/2;
        convmap=C(Nout+1:Nout+1024,Nout+1:Nout+1024);
        [Cl,l]=get_angular_spec(map,map,pixsize);
        [Clc]=get_angular_spec(convmap,convmap,pixsize);
        bl_arr(isim,:)=Clc./Cl;
    end    
    bldat.pixsize=pixsize;
    bldat.bl(ifield).bl=squeeze(mean(bl_arr));
    bldat.bl(ifield).l=l;
    bldat.bl(ifield).blerr=squeeze(std(bl_arr));
    
    bl_arr=zeros(nsim,29);
    for isim=1:nsim
        disp(sprintf('ifield%d,isim%d',ifield,isim));
        map=randn(1024);map=map-mean(map(:));
        C=conv2(map,bigpsf);
        Nout=(size(C,1)-size(map,1))/2;
        convmap=C(Nout+1:Nout+1024,Nout+1:Nout+1024);
        [Cl,l]=get_angular_spec(map,map,7);
        [Clc]=get_angular_spec(convmap,convmap,7);
        bl_arr(isim,:)=Clc./Cl;
    end    
    bldat.bl_big(ifield).bl=squeeze(mean(bl_arr));
    bldat.bl_big(ifield).l=l;
    bldat.bl_big(ifield).blerr=squeeze(std(bl_arr));
    
    save(strcat(savedir,'bldat'),'bldat');
end
return