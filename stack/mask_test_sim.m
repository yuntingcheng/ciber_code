flight=40030;
inst=1;
pixsize=7;
savedir='/Users/ytcheng/ciber/doc/20170617_Stacking/plots/cats/';
%%
psfsig=1;

mst_arr=10:1:17;
Ist_arr=10.^(-mst_arr/2.5);
Nst_arr=(Ist_arr.^-1);
Nst_arr=Nst_arr./Nst_arr(1);

strmap=zeros(1024);
for i=numel(mst_arr):-1:2
    mst=mst_arr(i);
    Ist=Ist_arr(i);
    Nst=round(Nst_arr(i));
    pos_arr=unifrnd(1,1024,[Nst,2]);
    posdat(i).mst=mst;
    posdat(i).Ist=Ist;
    posdat(i).Nst=Nst;
    posdat(i).pos_arr=pos_arr;
    for j=1:Nst
        radmap=make_radius_map(randn(1024),pos_arr(j,1),pos_arr(j,2));
        psf_arr=Ist./sqrt(2*pi*psfsig^2)*exp(-radmap./2./psfsig^2);
        strmap=strmap+psf_arr;
    end
    [Cl,l,~,~,~,~,Cl2d]=get_angular_spec(strmap,strmap,pixsize);
    loglog(l,l.*(l+1).*Cl./2./pi,'color',get_color(i-1));hold on
    drawnow
end
%%%
alpha=-15;
beta=250;

mask=ones(1024);
for i=1:numel(mst_arr)-1
    mst=posdat(i).mst;
    Ist=posdat(i).Ist;
    Nst=posdat(i).Nst;
    pos_arr=posdat(i).pos_arr;
    rad=alpha.*mst+beta;
    rad=rad./7;
    
    if rad<1
        rad=1;
    end

    for j=1:Nst
        radmap=make_radius_map(randn(1024),pos_arr(j,1),pos_arr(j,2));
        sp = find (radmap < rad);
        mask(sp)=0;
    end
    
    mmap=strmap.*mask;
    mmap=mmap-mean(mmap(find(mask)));
    [Cl,l]=get_angular_spec(mmap,mmap,pixsize);
    loglog(l,l.*(l+1).*Cl./2./pi,'o--','color',get_color(i));hold on
    drawnow    
end

%%% get bl
for frac=[0,0.25,0.5]
radmap=make_radius_map(randn(1025),513+frac,513+frac);
psf_arr=1./sqrt(2*pi*psfsig^2)*exp(-radmap./2./psfsig^2);
[bl,l]=get_angular_spec(psf_arr,psf_arr,pixsize);

loglog(l,l.*(l+1).*bl.*1e-11./2./pi,'-','color','k','LineWidth',3);hold on
end
%%
xlim([1e2,2e5]);
xlabel('$\ell$','interpreter','latex','fontsize',18)
ylabel('$\ell(\ell+1)C_\ell/2\pi(nW^2 m^{-4} sr^{-2})$',...
        'interpreter','latex','fontsize',18)
savename=strcat(savedir,'mask_test_sim');
print(savename,'-dpng');