function A = get_profile_mkk(flight,inst,ifield,dx,binedges)

mypaths=get_paths(flight);
psfdir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(psfdir,'fitpsfdat'),'fitpsfdat');

bestparam = fitpsfdat(ifield).bestparam;
beta = bestparam(1);
rc = bestparam(2);

radmap = make_radius_map(zeros(2*dx+1),dx,dx).*0.7;
psfmap = (1 + (radmap/rc).^2).^(-3.*beta./2);
psfmap = psfmap./sum(psfmap(:));
psfmap_small = psfmap(dx+1-50:dx+1+50,dx+1-50:dx+1+50);

[xx,yy] = meshgrid(-20:20,-20:20);
psf_pix = (10 - abs(xx)).*(10 - abs(yy));
psf_pix(abs(xx)>=10 | abs(yy) >= 10) = 0;
psf_pix = psf_pix ./ sum(psf_pix(:));

nbins = numel(binedges) - 1;
A = zeros(nbins);
for ibin = 1:nbins
    radmap = make_radius_map(zeros(2*dx+1),dx+1,dx+1);
    sigmap = zeros(2*dx+1);
    sigmap(radmap >= binedges(ibin) & radmap < binedges(ibin+1)) = 1;
    obsmap = conv2(sigmap, psfmap_small,'same');
    obsmap = conv2(obsmap, psf_pix,'same');
    profile = radial_prof(obsmap,ones(2*dx+1),dx+1,dx+1,1,0,'binedges',binedges);
    A(:,ibin) = profile.prof;
end
return