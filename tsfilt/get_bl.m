function get_bl(flight,inst)
mypaths=get_paths(flight);
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

nsim=50;
pixsize=0.7;
dx = 30;

radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*pixsize;
for ifield=4:8
    beta = fitpsfdat(ifield).psfmodel.beta_best;
    rc = fitpsfdat(ifield).psfmodel.rc_best;
    norm = fitpsfdat(ifield).psfmodel.norm;
    psfmap =norm.*(1 + (radmap/rc).^2).^(-3.*beta./2);
    bl_arr=zeros(nsim,29);
    for isim=1:nsim
        fprintf('ifield%d,isim%d\n',ifield,isim);
        map=randn(10240);
        C=conv2(map,psfmap);
        Nout=(size(C,1)-size(map,1))/2;
        convmap=C(Nout+1:Nout+10240,Nout+1:Nout+10240);
        map = rebin_map_coarse(map,10);
        convmap = rebin_map_coarse(convmap,10);
        [Cl,l]=get_angular_spec(map,map,pixsize*10);
        [Clc]=get_angular_spec(convmap,convmap,pixsize*10);
        bl_arr(isim,:)=Clc./Cl;
    end    
    bldat.pixsize=pixsize*10;
    bldat.bl(ifield).bl=squeeze(mean(bl_arr));
    bldat.bl(ifield).l=l;
    bldat.bl(ifield).blerr=squeeze(std(bl_arr));
    
    save(strcat(savedir,'bldat'),'bldat');
end
return