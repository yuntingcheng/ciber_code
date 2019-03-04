pixsize = 0.7;
dx = 1200;
beta = 1;
rc = 2.5;
radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1).*0.7;
psfmap_true = (1 + (radmap./rc).^2).^(-3.*beta./2);    
profile = radial_prof(psfmap_true,ones(2*dx+1),dx+1,dx+1,1,32);
ptrue_arr=(profile.prof)./profile.prof(1);
r_arr = profile.r.*pixsize;
%%
Nsims = 100;
src_coord=dx+0.5+10*rand(2,Nsims);
    stamp=zeros(2*dx+1,2*dx+1);
for isim = 1:Nsims
    xsrc=src_coord(1,isim);
    ysrc=src_coord(2,isim);
    radmap = make_radius_map(zeros(2*dx+10),xsrc,ysrc).*0.7;
    psfmap = (1 + (radmap./rc).^2).^(-3.*beta./2);    
    psfmap_coarse=rebin_map_coarse(psfmap,10);
    psfmap_fine=imresize(psfmap_coarse,10,'method','nearest');
    stamp=stamp+psfmap_fine(round(xsrc)-dx:round(xsrc)+dx,...
                            round(ysrc)-dx:round(ysrc)+dx);
end
stack=stamp./Nsims;

profile = radial_prof(stack,ones(2*dx+1),dx+1,dx+1,1,32);
pstack_arr=(profile.prof)./profile.prof(1);
%%
tophat = zeros(size(stack));
tophat(1196:1206, 1196:1206) = 1;
[~,psfr] = deconvblind(stack,tophat);

profile = radial_prof(psfr,ones(2*dx+1),dx+1,dx+1,1,32);
pr_arr=(profile.prof)./profile.prof(1);

profile = radial_prof(J,ones(2*dx+1),dx+1,dx+1,1,32);
pj_arr=(profile.prof)./profile.prof(1);
%%
loglog(r_arr,ptrue_arr);hold on
loglog(r_arr,pstack_arr);
loglog(r_arr,pr_arr);
loglog(r_arr,pj_arr);
ylim([1e-3,1.1])
legend({'true','stack','recovered','J'})