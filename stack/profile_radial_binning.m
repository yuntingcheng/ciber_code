function [prof15,prof100] = profile_radial_binning(prof25,w,sp100)
    % binning the raidial bins such that the cov is invertible
    % also return the avg profile > 100 arcsec
    
    prof15 = zeros([1,15]);
    prof15(2:end-1) = prof25(7:19);
    prof15(1) = sum(prof25(1:6).*w(1:6))./sum(w(1:6));
    prof15(end) = sum(prof25(20:25).*w(20:25))./sum(w(20:25));
    
    prof100 = sum(prof25(sp100).*w(sp100))./sum(w(sp100));
return
[]