function [sigmap, shotmap]=sigmap_from_mz14(npix,pixscale)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Generate Gaussian random signal map and shot noise map from MZ PS.
%Input:
%   -npix: # of pixels (1024)
%   -pixscale: pixsize in arcsec ( 7'')
%Output:
%   -sigmap
%   -shotmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~,l,~,Clsig,Clshot]=sigCl_extrap_mz14;
    l2d = get_l(npix,npix,pixscale);
    
    mapfft2 = interp1(l,sqrt(Clsig),l2d,'spline');
    sigmap = real(map_from_FFT2(mapfft2,pixscale));
    sigmap = sigmap - mean(sigmap(:));

    mapfft2 = interp1(l,sqrt(Clshot),l2d,'spline');
    shotmap = real(map_from_FFT2(mapfft2,pixscale));
    shotmap = shotmap - mean(shotmap(:));
return
