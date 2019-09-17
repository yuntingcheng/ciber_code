function [prof_arr,profpsf_arr,r_arr] =  PSF_stacked_profile(flight,inst,ifield)
mypaths=get_paths(flight);
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');
beta=fitpsfdat(ifield).psfmodel.beta_best;
rc=fitpsfdat(ifield).psfmodel.rc_best;
norm=fitpsfdat(ifield).psfmodel.norm;

dx=1200;
Nlarge =dx+9;
radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
Imap_large = norm .* (1 + (radmap/rc).^2).^(-3.*beta./2);

stack = zeros(2*dx+1);
for i=1:10
    for j=1:10
        stamp = Imap_large(i:i+2*dx+9,j:j+2*dx+9);
        stamp = rebin_map_coarse(stamp,10);
        stamp = stamp.*100;
        stamp = imresize(stamp,10,'nearest');
        stamp = stamp(1211-i-dx:1211-i+dx,1211-j-dx:1211-j+dx);
        stack = stack + stamp;
    end
end

nbins = 25;
profile = radial_prof(stack,ones(2*dx+1),dx+1,dx+1,1,nbins);
prof_arr = profile.prof./profile.prof(1);
rbinedges = profile.binedges;
r_arr = binedges2bins(rbinedges).*0.7;

profile = radial_prof(Imap_large(Nlarge+1-dx:Nlarge+1+dx,...
    Nlarge+1-dx:Nlarge+1+dx),ones(2*dx+1),dx+1,dx+1,1,nbins);
profpsf_arr = profile.prof./profile.prof(1);
return